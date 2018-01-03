import sh
from corral import cli, conf, db
from . import models

class SqliteBrowser(cli.BaseCommand):
    """Opens sqlitebrowser GUI to look SQL DB entries
    """
    def handle(self):
        if conf.settings.CONNECTION.startswith("sqlite:///"):
            path = conf.settings.CONNECTION.replace("sqlite:///", "")
            sh.sqlitebrowser(path)

#~ class SetupImagesCount(cli.BaseCommand):
    #~ """Populates DB with first static entries
    #~ """
    #~ def handle(self):
        #~ image = models.Images()
        #~ image.path = '/home/bruno/Data/RBpipelin/images/img-1'

        #~ with db.session_scope() as session:
            #~ session.add(image)




