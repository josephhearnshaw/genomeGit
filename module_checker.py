import pip

#Required modules
pakgs = ['pytabix', 'pyfaidx']
for package in pakgs:
    try:
        import package
    except ImportError, e:
        pip.main(['install', package])
