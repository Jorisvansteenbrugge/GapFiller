from setuptools import setup, find_packages

setup(name='Hydraslayer',
      version="1.2.0",
      description=u"Program that can 'slay' gaps in dna sequences using (preferably) long reads",
      author=u'Joris van Steenbrugge',
      author_email='joris.vansteenbrugge@wur.nl',
      license="MIT",
      packages=['Hydraslayer'],
      scripts=["scripts/Hydraslayer"]
      )
