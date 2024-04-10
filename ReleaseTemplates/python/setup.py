from setuptools import setup
from setuptools.dist import Distribution

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False
    
setup(
    name='neurannparser',
    version='0.0.1',    
    description='Parser for NEURANN Chromosomes and Networks',
    url='https://github.com/NEURANN/NetworkParser',
    author='AyOhEe',
    author_email='amacd1124@gmail.com',
    license='MIT',
    packages=['neurannparser'],
    package_data= {'neurannparser' : ['NetworkParser.dll']},
    distclass=BinaryDistribution,

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research'
        'License :: OSI Approved :: MIT License',  
        'Operating System :: Microsoft :: Windows',        
        'Programming Language :: Python :: 3.9',        
        'Programming Language :: Python :: 3.10',        
        'Programming Language :: Python :: 3.11',        
    ],
)