import setuptools

with open('readme.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyvivc',                     
    version='1.3',                       
    author='Stephen Checkley',                    
    author_email='scheckley@gmail.com',
    url='https://github.com/scheckley/pyvivc',
    description="numerical deconvolution and convolution methods working for inequal and incompatible timepoints between impulse and response curves for application in IVIVC level A",
    long_description=long_description,      
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),    
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],                                      
    python_requires='>=3.9',      
    py_modules=['pyvivc','newrange'],            
    package_dir={'':'pyvivc/src'},   
    install_requires=['scipy', 'pandas', 'numpy'],
    project_urls={
        'Bug Reports': 'https://github.com/scheckley/pyvivc'
    }   
)
