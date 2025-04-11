from setuptools import setup, find_packages

setup(
    name="deepsmirud",
    version="0.1.2",
    keywords=["pip", "deepsmirud"],
    description="deepsmirud",
    long_description="predicting types of small molecule regulation on microRNA expression",
    license="GPL-3.0",

    url="https://github.com/2003100127/deepsmirud",
    author="Jianfeng Sun",
    author_email="jianfeng.sunmt@gmail.com",

    packages=find_packages(),
    include_package_data=True,
    platforms="any",
    python_requires='>=3.11,<3.12',
    install_requires=[
        'pandas',
        'biopython',
        'pyfiglet',
        'click',
        'rdkit',
        'numpy==1.24.3',
        'tensorflow==2.14',
    ],
    entry_points={
        'console_scripts': [
            'deepsmirud=deepsmirud.predict:run',
            'deepsmirud_download=deepsmirud.predict:download',
        ],
    }
)