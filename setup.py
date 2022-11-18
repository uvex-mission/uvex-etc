from setuptools import setup, find_packages

setup(
    name='uvex_etc',
    version='0.1dev',
    packages=find_packages(),
    install_requires=[
        'requests',
        'importlib; python_version == ">3.8"',
    ],
)
