
from setuptools import setup

if __name__ == "__main__":
    setup(
            name='chemflow',
            url='https://github.com/turnerluke/chemflow',
            author='Turner LUke',
            author_email='tluke@wisc.edu',
            # Needed to actually package something
            packages=['measure'],
            # Needed for dependencies
            install_requires=['numpy', 'scipy', 'pint'],
            # *strongly* suggested for sharing
            version='0.1',
            # The license can be anything you like
            license='MIT',
            description='chemflow a Python package for quick sequential chemical process design calculations.'
            # long_description=open('README.txt').read(),
    )
