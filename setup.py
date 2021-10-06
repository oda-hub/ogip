from setuptools import setup
import sys

setup_requires = ['setuptools >= 30.3.3']

if {'pytest', 'test', 'ptr'}.intersection(sys.argv):
    setup_requires.append('pytest-runner')


setup(description="ogip",
      long_description=open('README.md').read(),
      version='0.1.3',
      include_package_data=True,
      setup_requires=setup_requires)
