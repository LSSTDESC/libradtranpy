from setuptools import setup

setup(
   name='libradtranpy',
   version='0.2.0',
   author='Sylvie Dagoret-Campagne',
   author_email='sylvie.dagoret-campagne@ijclab.in2p3.fr',
   packages=['libradtranpy', 'test'],
   scripts=['bin/script.sh'],
   url='https://github.com/LSSTDESC/libradtranpy',
   license='LICENSE',
   description='python wrapper for libradtran',
   long_description=open('README.md').read(),
   install_requires=[
       "Django >= 1.1.1",
       "pytest",
   ],
)

