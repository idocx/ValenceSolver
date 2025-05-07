from setuptools import setup, find_packages

__author__ = 'Tanjin He'
__maintainer__ = 'Tanjin He'
__email__ = 'tanjin_he@berkeley.edu'

if __name__ == "__main__":
    setup(name='ValenceSolver',
          version=1.0,
          author="Tanjin He",
          author_email="tanjin_he@berkeley.edu",
          license="MIT License",
          packages=find_packages(),
          package_data={'ValenceSolver': ['core/analysis/icsd_bv.yaml']},
          install_requires=[
              "sympy",
              "unidecode",
              "pymatgen",
              "pulp",
          ],
          
          zip_safe=False)

