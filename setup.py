import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyvivc',                     
    version='1.0',                       
    author='Stephen Checkley',                    
    author_email='scheckley@gmail.com',
    url='',
    description="numerical deconvolution and convolution methods working for inequal and incompatible timepoints between impulse and response curves for application in IVIVC level A",
    long_description=long_description,      
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),    
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: GPL :: Level 3',
        'Operating System :: OS Independent',
    ],                                      
    python_requires='>=3.9',      
    py_modules=['pyvivc'],            
    package_dir={'':'pyvivc/src'},   
    install_requires=['scipy, pandas, numpy'],
    project_urls={
        'Bug Reports': ''
    }   
)