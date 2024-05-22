from setuptools import setup,find_packages

setup(name="scattersim",
      description="Meso scale SAXS Simulations",
      version='1.0',
      author='Kevin Yager',
      license='MIT',
      python_requires='>=3.8',
      install_requires=['numpy',
                        'scipy', 
                        'matplotlib', 
                        ],
      extras_require = {},
      packages=find_packages(),
      long_description=open('README.md').read(),
      long_description_content_type="text/markdown",
      classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows"
      ],
)