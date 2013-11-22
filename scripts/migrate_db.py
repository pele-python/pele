from pele.storage import database
import sqlalchemy
import sys

def from_0_to_1(connection, schema):
    ''' migrating from version 0 to 1
    
        fields for log product of frequencies and point group order were added
        to Minimum and TransitionState
    '''
    assert schema == 0
    print "migrating from database version 0 to 1"
    connection.execute("ALTER TABLE tbl_minima ADD fvib FLOAT;")
    connection.execute("ALTER TABLE tbl_minima ADD pgorder INTEGER;") 
    connection.execute("ALTER TABLE tbl_transition_states ADD fvib FLOAT;")
    connection.execute("ALTER TABLE tbl_transition_states ADD pgorder INTEGER;") 
    connection.execute("PRAGMA user_version = 1;")
    return 1

def from_1_to_2(connection, schema):
    assert schema == 1
    print "migrating from database version 1 to 2"
    connection.execute("""CREATE TABLE tbl_system_property (
    _id INTEGER NOT NULL, 
    property_name VARCHAR, 
    int_value INTEGER, 
    float_value FLOAT, 
    string_value VARCHAR, 
    pickle_value BLOB, 
    PRIMARY KEY (_id)
    );
    """)
    connection.execute("PRAGMA user_version = 2;")
    return 2

migrate_script = {0:from_0_to_1,
                  1:from_1_to_2
                  }
    
def migrate(db):
    engine = sqlalchemy.create_engine("sqlite:///%s"%db)
    res = engine.execute("PRAGMA user_version;")
    schema = res.fetchone()[0]
    res.close()
    
    print "current version:",schema
    print "newest version:",database._schema_version
    
    connection = engine.connect()
    while schema < database._schema_version:
        trans = connection.begin()
        try:
            schema = migrate_script[schema](connection, schema)
        except BaseException:
            print "failed to migrate database, rolling back changes"
            trans.rollback()
            raise
        trans.commit()
    connection.close()
    print "database is at newest version"
    
if __name__ == '__main__':
    if len(sys.argv) < 2 or "--help" in sys.argv or "-h" in sys.argv:
        print "usage:\npython migrate_db.py <database1> [<database2> ...]"
        print ""
        print "update pele database to the newest version"
        sys.exit()
    for dbfile in sys.argv[1:]:
        print dbfile
        migrate(dbfile)