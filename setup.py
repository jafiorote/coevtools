from setuptools import setup

setup(name='Coevtools',
      version='0.1',
      description='',
      url='https://github.com/jafiorote/coevtools',
      author='José Antonio Fiorote',
      author_email='jafiorote@gmail.com',
      license='GPL 3.0',
      packages=['coevtools'],
      install_requires=[
            'biopython',
            'python-Levenshtein',
            'numpy', 'scipy', 'Bio'],
      zip_safe=False)
