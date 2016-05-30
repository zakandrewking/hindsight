try:
    from setuptools import setup, Command
except:
    from distutils.core import setup, Command

setup(
    name='hindsight',
    author='Zachary King',
    packages=['hindsight'],
)
