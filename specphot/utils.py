import numpy as np
import astropy.units as u

class Spectrum:
    def __init__(self, wavelength, flux, flux_err, mask, wavelength_units='um', flux_units='erg/s/cm**2/AA'):
        """
        This class creates an object to store a spectrum. The spectrum is 
        automatically converted to Angstroms and f_lambda units in cgs.

        Parameters
        ----------
        wavelength : array-like
            The wavelength of the spectrum.
        flux : array-like
            The flux of the spectrum.
        flux_err : array-like
            The error of the spectrum.
        mask : array-like
            An array of 1s and 0s to mask out bad data points, where 1 denotes 
            a masked value.
        wavelength_units : str
            The units of the wavelength. Default is 'um'.
        flux_units : str
            The units of the flux. Default is 'erg/s/cm**2/AA'.
        """
        self.wavelength = (wavelength * u.Unit(wavelength_units)).to('AA').value
        self.flux = (flux * u.Unit(flux_units)).to('erg/s/cm**2/AA', equivalencies=u.spectral_density(wavelength*u.Unit(wavelength_units))).value
        if flux_err is not None:
            self.flux_err = (flux_err * u.Unit(flux_units)).to('erg/s/cm**2/AA', equivalencies=u.spectral_density(wavelength*u.Unit(wavelength_units))).value
        else:
            self.flux_err = np.zeros_like(self.flux)
        self.mask = mask
        self.wavelength_units = wavelength_units
        self.flux_units = flux_units


    def clean_spectrum(self, filter_min, filter_max):
        """
        This function cleans the spectrum by replacing nan values and/or masked 
        values with the median value. It also extrapolates the spectrum to the
        ends of the filter if the filter extends beyond the ends of the spectrum.
        """
        # Replace nan values with the median value
        clean_wave = self.wavelength.copy() 
        clean_flux = self.flux.copy()
        clean_error = self.flux_err.copy()
        
        bad_points = np.isnan(clean_flux) | (self.mask == 1)

        clean_flux[bad_points] = np.nanmedian(clean_flux)
        clean_error[bad_points] = np.nanmedian(clean_error)

        # Extrapolate the spectrum to the ends of the filter with the same
        # wavelength sampling as the spectrum
        if filter_min < self.wavelength[0]:
            extrap_wl = np.linspace(filter_min, self.wavelength[0], np.diff(self.wavelength)[0])
            extrap_flux = np.ones(len(extrap_wl)) * np.nanmedian(clean_flux)
            extrap_error = np.ones(len(extrap_wl)) * np.nanmedian(clean_error)
            clean_wave = np.concatenate((extrap_wl, clean_wave))
            clean_flux = np.concatenate((extrap_flux, clean_flux))
            clean_error = np.concatenate((extrap_error, clean_error))

        if filter_max > self.wavelength[-1]:
            extrap_wl = np.linspace(self.wavelength[-1], filter_max, np.diff(self.wavelength)[-1])
            extrap_flux = np.ones(len(extrap_wl)) * np.nanmedian(clean_flux)
            extrap_error = np.ones(len(extrap_wl)) * np.nanmedian(clean_error)
            clean_wave = np.concatenate((clean_wave, extrap_wl))
            clean_flux = np.concatenate((clean_flux, extrap_flux))
            clean_error = np.concatenate((clean_error, extrap_error))

        return clean_wave, clean_flux, clean_error

    def pass_through_filter(self, filter_wavelength, filter_transmission, filter_wavelength_units='um', photometry_units='erg/s/cm**2/AA'):
        """
        This function passes the spectrum through a filter curve in order to 
        compute the synthetic photometry. All NaN values are set to the median
        flux density and error of the spectrum. If the filter extends beyond the
        ends of the spectrum, the spectrum will be extrapolated to the ends of
        the filter with the median flux density and error.

        Parameters
        ----------
        filter_wavelength : array-like
            The wavelength of the filter curve.
        filter_transmission : array-like
            The transmission of the filter curve.
        filter_wavelength_units : str
            The units of the filter wavelength. Default is 'um'.
        photometry_units : str
            The units of the photometry. Default is 'erg/s/cm**2/AA'.

        Returns
        -------
        f_filter : float
            The synthetic photometry.
        e_filter : float
            The error on the synthetic photometry.
        """
        # Interpolate the filter transmission to the spectrum wavelength array
        filter_wavelength = (filter_wavelength * u.Unit(filter_wavelength_units)).to('AA').value
        filter_min = np.amin(filter_wavelength); filter_max = np.amax(filter_wavelength)
        transmission = np.interp(self.wavelength, filter_wavelength, filter_transmission, left=0.0, right=0.0)
        lam_center = np.average(filter_wavelength, weights=filter_transmission)

        # Clean up nan values and extrapolate the spectrum to the ends of the filter.
        clean_wave, clean_flux, clean_error = self.clean_spectrum(filter_min, filter_max)

        dl = np.concatenate((np.diff(clean_wave), np.array([clean_wave[-1] - clean_wave[-2]])))
        flam_filter = np.sum(clean_flux * transmission * dl) / np.sum(transmission * dl)
        
        # Compute the error
        flam_filter_error = np.sqrt(np.sum((clean_error * transmission)**2 * dl) / np.sum(transmission * dl))

        f_filter = (flam_filter * u.Unit('erg/s/cm**2/AA')).to(u.Unit(photometry_units), equivalencies=u.spectral_density(lam_center*u.AA)).value
        e_filter = (flam_filter_error * u.Unit('erg/s/cm**2/AA')).to(u.Unit(photometry_units), equivalencies=u.spectral_density(lam_center*u.AA)).value
        lam_center = (lam_center * u.AA).to(u.Unit(filter_wavelength_units)).value
        
        return f_filter, e_filter, lam_center

    