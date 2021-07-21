from setuptools import setup

setup(name='Docking Score Module',
      version='0.1',
      description='',
      url='https://github.com/jafiorote/dsm',
      author='José Antonio Fiorote',
      author_email='jafiorote@gmail.com',
      license='MIT',
      packages=['dsm'],
      install_requires=[
            'biopython',
            'python-Levenshtein',
            'numpy', 'scipy', 'Bio'],
      zip_safe=False)
