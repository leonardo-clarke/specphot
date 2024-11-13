from setuptools import setup, find_packages

setup(
    name='specphot',
    version='0.1.0',
    description='This Python package provides tools for comparison of spectra and photometry.',
    author='Leonardo Clarke',
    author_email='leoclarke@astro.ucla.edu',
    url='https://github.com/leonardo-clarke/specphot',
    packages=find_packages(),
    install_requires=[],  # list dependencies here or load from requirements.txt
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
)