import pip

# Required modules
pakgs = ['pytabix==0.1', 'pyfaidx==0.5.5.2']
for package in pakgs:
    try:
        import package
    except ImportError, e:
        pip.main(['install', package])
