from setuptools import setup, Command, find_packages

setup(
    name='hindsight',
    author='Zachary King',
    python_requires='>=2.7,<3.0',
    packages=find_packages(),
    install_requires=[
        'cobra>=0.5,<0.6',
        'jupyter',
        'numpy',
        'scipy',
        'pandas',
        'matplotlib_venn',
        'seaborn',
        'sklearn',
    ],
)
