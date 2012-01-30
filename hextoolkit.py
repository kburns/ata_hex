"""
Commonly used hex analysis code.

Author: Keaton J. Burns <keaton.burns@gmail.com>
Affiliation: UC Berkeley

"""


import os
import os.path
from os.path import join
import numpy as np
import sqlite3


# Set things up so we can dump numpy array values right into sqlite databases.
sqlite3.register_adapter(np.int32, int)
sqlite3.register_adapter(np.float32, float)
sqlite3.register_adapter(np.string_, str)


def getdbpath():
    """Retrieve path to database file"""
    if 'SQUINTDBPATH' in os.environ:
        return os.environ['SQUINTDBPATH']
    return '/ataarchive/scratch/hexproc/squint.db'


def infopath(*args):
    return join(os.path.dirname(__file__), *args)


def antnum(astr):
    """Return ATA antenna number given antenna name"""
    telnum = np.genfromtxt(infopath('telnum.txt'), dtype='|S2', comments='#')
    return telnum.tostring().index(astr)/2 + 1


def antname(anum):
    """Return ATA antenna name given antenna number"""
    telnum = np.genfromtxt(infopath('telnum.txt'), dtype='|S2', comments='#')
    return telnum[anum - 1]


def sqlitedb_to_ndarray(fname, table='data', str_length=20):
    """Load table from a table in a sqlite3 .db file to a numpy ndarray"""
    
    connection = sqlite3.connect(fname)
    cursor = connection.cursor()
    cursor.execute('select * from ' + table)
    
    data = sqlitecursor_to_ndarray(cursor, str_length=str_length)
    
    # closing the connection
    connection.close()
    
    return data
	
	
def sqlitecursor_to_ndarray(cursor, str_length = 20):
    """Load table from sqlite3 cursor selection to a numpy ndarray"""

    # Get data types
    types = []
    data = cursor.fetchall()
    if data == []: return None
    for i in xrange(np.size(data[0])):
        if type(data[0][i]) == unicode: 
            # Take longest string over all rows
            istrlength = str_length
            for j in data:
                if len(j[i]) > istrlength: istrlength = len(j[i])
            types.append('S%s' % istrlength)
        if type(data[0][i]) == float: types.append('float')
        if type(data[0][i]) == int: types.append('int')

    # Get column names
    varnm = [i[0] for i in cursor.description]

    # Autodetected dtype
    dtype = zip(varnm, types)
    data = np.array(data, dtype=dtype)

    return data
    
    
def atatojday(atadate):
    """
    atatojday
    =========
    
    PURPOSE:
        Converts data-databounds date strings into Julian days
        
    CALLING SEQUENCE:
        atatojday(datestr)
      
        Sample atadate: '10Sep08:11:25:55.0'
      
    """ 
    
    # Split off semicolons and pull out components
    atadate = atadate.split(':')
    
    date = atadate[0]
    hour = float(atadate[1])
    min = float(atadate[2])
    sec = float(atadate[3])
    
    year = int('20' + date[0:2])
    month = date[2:-2]
    day = int(date[-2:])
    
    if month == 'Jan': month = 1
    elif month == 'Feb': month = 2
    elif month == 'Mar': month = 3
    elif month == 'Apr': month = 4
    elif month == 'May': month = 5
    elif month == 'Jun': month = 6
    elif month == 'Jul': month = 7
    elif month == 'Aug': month = 8
    elif month == 'Sep': month = 9
    elif month == 'Oct': month = 10
    elif month == 'Nov': month = 11
    elif month == 'Dec': month = 12
    
    # Construct Julian Day (by Wikipedia's algorithm)
    a = (14 - month) / 12
    y = year + 4800 - a
    m = month + 12*a - 3
    
    JDN = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045
    JD = JDN + (hour - 12) / 24. + min / 1440. + sec / 86400.
    return JD
    

def feedID(antnum, julday):
    """Return feed number for antenna antum on day julday"""
  
    # Read in feedswapjd.txt file
    feedlog = np.genfromtxt(infopath('feedswapjd.txt'), dtype=('f' + ',f' * 43), delimiter=',')
    
    # Find all swaps before requested julday, take ant from latest one
    afterswitch = np.where(feedlog['f0'] < julday * np.ones(np.shape(feedlog['f0'])))
    jloc = np.where(feedlog['f0'] == np.max(feedlog['f0'][afterswitch]))
     
    return feedlog[jloc][0][antnum + 1]
   

def gaussread(path):
    """
    gaussread
    =========
    
    PURPOSE:
        Reads data textfiles including:
            data-gaussfits, data-sefd, data-sinfo, data-archsummpath, data-databounds
        Reduces to numpy ndarrays
    
    CALLING SEQUENCE:
        [gread, eread, info, squint] = gaussread(path)
    
    INPUTS:
        path    :=  path to folder containing hex reduction txts
    
    """
    
    
    ##################
    ## Read in Data ##
    ##################
    
    # Read in gaussian fit file, bail if empty
    GDTYPES = 'i,S1,i,i,f,f,f,f,f,f,f,f,f,f,f'
    GNAMES = ('ANT,POL,Npts,,XiSq,AMP,AMPuc,OffAz,OffAzuc,OffEl,'
              'OffEluc,WidthAz,WidthAzuc,WidthEl,WidthEluc')
    gread = hexfromtxt(join(path, 'data-gaussfits.txt'), dtype=GDTYPES, names=GNAMES, colnum=15)
    if gread == None: return 'GAUSSFITS_READ_FAILURE', 0, 0, 0
    
    # Read in SEFD file, intend to pass zeros if empty
    EDTYPES = 'i,S1,f,f,f,f,f'
    ENAMES = 'Ant,Pol,Avg-Amp,Amp-RMS,Avg-Pha,Pha-RMS,SEFD'
    eread = hexfromtxt(join(path, 'data-sefd.txt'), dtype=EDTYPES, names=ENAMES, skip_header=2)
    if eread == None:
        eread = 'SEFD_READ_FAILURE'
        print '   No SEFD data found.'

    # Read in run information
    sinfofile = open(join (path, 'data-sinfo.txt'), 'r')
    sinforead = sinfofile.readlines()
    sinfofile.close()
    
    sinforead = [i.strip().split() for i in sinforead]
    source = sinforead[0][1]
    freq = sinforead[1][1]
    flux = sinforead[2][1]
    
    # Read in path to archive summ path, will use as unique run identifier
    archfile = open(join (path, 'data-archsummpath.txt'), 'r')
    archread = archfile.readlines()
    archfile.close()
    archsummpath = archread[0].strip()
    
    # Read in start and end times, convert to Julian dates
    boundsfile = open(join (path, 'data-databounds.txt'), 'r')
    boundsread = boundsfile.readlines()
    boundsfile.close()

    boundsread = [i.strip().split() for i in boundsread]
    tstart = atatojday(boundsread[0][1])
    #tend = atatojday(boundsread[1][1])

    info = (tstart, source, freq, flux, archsummpath)
        
    
    ###########################
    ## Reduce to useful data ##
    ###########################
    
    # Find out how many x,y center pairs there are
    sqpairs = 0
    for i in xrange(1,43):
        if np.size(np.where(gread['ANT'] == i)) == 2: sqpairs += 1
    
    if sqpairs == 0:
        print '   No squint pairs found.'
        # We'll still add the run to the DB.
        #return 'NO_SQUINT_PAIRS', 0, 0, 0
    
    # Define types and names for the squint ndarray
    SDTYPES = 'i S2 f f f f f f f f f f f f f f f'.split ()
    SNAMES = ['antnum', 'antname', 'feed', 'sefd', 'sumchisq',
              'squintaz', 'squintaz_uc', 'squintel', 'squintel_uc',
              'x_width_az', 'x_width_az_uc', 'x_width_el', 'x_width_el_uc',
              'y_width_az', 'y_width_az_uc', 'y_width_el', 'y_width_el_uc']
    squint = np.zeros(sqpairs, dtype = zip(SNAMES, SDTYPES))
    
    # Calculate squints
    j = 0
    for i in xrange(1,43):
        antloc = np.where(gread['ANT'] == i)
        
        # Pick antennas with x and y gaussfit data
        if np.size(antloc) == 2:
            squint[j]['antnum'] = i
            squint[j]['antname'] = antname(i)
            
            # Add feed info
            squint[j]['feed'] = feedID(i, tstart)
            
            # Squint defined from x to y, in arcminutes
            xloc = np.where(gread[antloc]['POL'] == 'x')
            yloc = np.where(gread[antloc]['POL'] == 'y')
            squint[j]['squintaz'] = (gread[antloc][yloc]['OffAz'] - gread[antloc][xloc]['OffAz']) / 60.0
            squint[j]['squintaz_uc'] = np.sqrt (gread[antloc][yloc]['OffAzuc']**2 +
                                                gread[antloc][xloc]['OffAzuc']**2) / 60.0
            squint[j]['squintel'] = (gread[antloc][yloc]['OffEl'] - gread[antloc][xloc]['OffEl']) / 60.0
            squint[j]['squintel_uc'] = np.sqrt (gread[antloc][yloc]['OffEluc']**2 +
                                                gread[antloc][xloc]['OffEluc']**2) / 60.0

            squint[j]['x_width_az'] = gread[antloc][xloc]['WidthAz']
            squint[j]['x_width_az_uc'] = gread[antloc][xloc]['WidthAzuc']
            squint[j]['x_width_el'] = gread[antloc][xloc]['WidthEl']
            squint[j]['x_width_el_uc'] = gread[antloc][xloc]['WidthEluc']
            
            squint[j]['y_width_az'] = gread[antloc][yloc]['WidthAz']
            squint[j]['y_width_az_uc'] = gread[antloc][yloc]['WidthAzuc']
            squint[j]['y_width_el'] = gread[antloc][yloc]['WidthEl']
            squint[j]['y_width_el_uc'] = gread[antloc][yloc]['WidthEluc']
            
            # Take SEFD as geometric mean of x and y for each antenna, if present
            if eread == 'SEFD_READ_FAILURE':
                squint[j]['sefd'] = 0.0
            else:
                eloc = np.where(eread['Ant'] == i)
                if np.size(eloc) != 0:
                    squint[j]['sefd'] = np.exp(np.log(eread[eloc]['SEFD']).mean())
                else: 
                    squint[j]['sefd'] = 0.0

            squint[j]['sumchisq'] = gread[antloc][xloc]['XiSq'] + gread[antloc][yloc]['XiSq']
            
            j += 1
            
    return gread, eread, info, squint
    
    
def gausstosql(path, replacedups=True):
    """
    gausstosql
    =========
    
    PURPOSE:
        Adds squint array data from gaussread to the squint SQL database.
    
    CALLING SEQUENCE:
        gausstosql(path, replacedups=True)
    
    INPUTS:
        path        :=  path to folder containing hex reduction txts
        replacedups :=  whether information for hex runs that have already
                        been entered into the database should replace the
                        preexisting information. If False, the new information
                        will be ignored.
    """
    
    # Run gaussread
    [gread, eread, info, squint] = gaussread(path)
    
    # Bail if fail
    if gread == 'GAUSSFITS_READ_FAILURE':
        print 'GAUSSTOSQL: Failed to read data-gaussfits.txt, closing...'
        return
    # If 'NO_SQUINT_PAIRS' flagging is activated in gaussread, this will catch
    if gread == 'NO_SQUINT_PAIRS':
        print 'GAUSSTOSQL: No squint pairs found, closing...'
        return
    
    # Connect to sqlite database
    connection = sqlite3.connect(getdbpath ())
    cursor = connection.cursor()
    
    # Check to see if already in database
    ASPtest = info[4]
    cursor.execute('SELECT rid FROM runs WHERE archsummpath = ?', 
                   (ASPtest, ))
    matches = cursor.fetchall ()

    if len (matches) > 0:
        rid = matches[0][0]

        if not replacedups:
            print 'GAUSSTOSQL: Run', ASPtest, 'already entered into database; dropping new data'
            connection.close ()
            return

        print 'GAUSSTOSQL: Run', ASPtest, 'already entered into database; replacing old data'
        cursor.execute('DELETE FROM runs WHERE rid = ?', (rid, ))
        cursor.execute('DELETE FROM obs WHERE rid = ?', (rid, ))

    # Insert data. NULL values will be replaced with automatic object ids.
    cursor.execute('INSERT INTO runs VALUES (NULL,?,?,?,?,?)', info)
    runID = cursor.lastrowid
    cursor.executemany ('INSERT INTO obs VALUES (NULL,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                        ((runID, ) + tuple(row) + (0,) for row in squint))
    
    # Finish up!
    connection.commit()
    connection.close()
    
    
def resetsql():
    """Reset sql database, backup old database to squint.db.old"""
 
    # Backup database
    dbpath = getdbpath()
    try:
        os.rename(dbpath, dbpath + '.old')
    except OSError, e:
        # Ignore the error if the database doesn't yet exist.
        if e.errno != 2:
            raise

    # Connect to sqlite database
    connection = sqlite3.connect(dbpath)
    cursor = connection.cursor()
    
    # Create table to hold observations
    sql_cmd = """CREATE TABLE obs (oid INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, 
                                   rid int,
                                   antnum int,
                                   antname text,
                                   feed float,
                                   sefd float,
                                   sumchisq float,
                                   squintaz float,
                                   squintaz_uc float,
                                   squintel float,
                                   squintel_uc float,
                                   x_width_az,
                                   x_width_az_uc,
                                   x_width_el,
                                   x_width_el_uc,
                                   y_width_az,
                                   y_width_az_uc,
                                   y_width_el,
                                   y_width_el_uc,
                                   flag int);"""
    cursor.execute(sql_cmd)
    
    # Create table to hold runs
    sql_cmd = """CREATE TABLE runs (rid INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, 
                                   date float,
                                   source text,
                                   freq float,
                                   flux float,
                                   archsummpath text);"""
    cursor.execute(sql_cmd) 
    
    # Close and exit
    connection.commit()
    connection.close()

    
def buildsql(rootdir):
    """
    buildsql
    ========
    
    PURPOSE:
        Adds to squint.db by sorting through all subdirectories of given root directory
        NOTE: Does not reset database, instead adds to current database
              To reset use resetsql()

    CALLING SEQUENCE:
        buildsql(rootdir)
    
    INPUTS:
        rootdir     :=  root directory to add from
    """
   
    # Find all folders with data-gaussfits.txt

    def load(unused, dirname, filenames):
        if 'data-gaussfits.txt' in filenames:
            print dirname
            gausstosql(dirname)

    os.path.walk(rootdir, load, None)
    
    
def flagsql(wherecmd='', add=True, num=1):
    """
    flagsql
    =======
    
    PURPOSE:
        Flags data in squint.db to exlcude from plotting
        Reset flags using    flagsql(num=0, add=False)
        
    CALLING SEQUENCE:
        flagsql(wherecmd='', add=True, num=1)
        
    INPUTS:
        wherecmd    :=  SQL WHERE command specifying data to be flagged
                            'WHERE ...'
        add         :=  Boolean, True: add num to current flag, False: set flag to num
        num         :=  flag number:
  
                            1   -   Squintaz mag
                            2   -   Squintaz_uc mag
                            4   -   Squintel mag
                            8   -   Squintel_uc mag
                            16  -   Sumchisq mag
                            32  -   SEFD mag
        
    """
    
    # Connect to sqlite database
    connection = sqlite3.connect(getdbpath())
    cursor = connection.cursor()
    
    if add:
        sql_cmd = 'UPDATE obs SET flag=flag+%i ' %num + wherecmd
    else:
        sql_cmd = 'UPDATE obs SET flag=%i ' %num + wherecmd
        
    cursor.execute(sql_cmd)
    
    connection.commit()
    connection.close()
    
    
def querysql(taglist, wherecmd='', exclude_flagged=True):
    """
    querysql
    ========
    
    PURPOSE:
        Query the SQL database, properly building required derived values
        
    CALLING SEQUENCE:
        data = querysql(taglist, wherecmd='', exclude_flagged=True)
        
    INPUTS:
        taglist     :=  list of tags to retrieve
        wherecmd    :=  'WHERE ...' command for specifying data in sql query
        exclude_flagged :=  remove datapoints which have been flagged

    TAG LIST:
        In SQL database:
            'date'
            'source'
            'freq'
            'flux'
            'archsummpath'      (not physical)
            'rid'               (not physical)
            'antnum'
            'antname'
            'feed'
            'feedrev'               Decimal feed revision
            'sefd'
            'sumchisq'
            'squintaz'
            'squintaz_uc'
            'squintel'
            'squintel_uc'
            'x_width_az'
            'x_width_az_uc'
            'x_width_el'
            'x_width_el_uc'
            'y_width_az'
            'y_width_az_uc'
            'y_width_el'
            'y_width_el_uc'
            'flag'
            
        Derived from SQL database:
            'antfeed'
            'antfeedrev'
            'squintmag'
            'squintmag_uc'
            'squintangle'
            'squintangle_uc'
            
    """
    
    
    # Connect to sqlite database
    connection = sqlite3.connect(getdbpath())
    cursor = connection.cursor()
    
    # Query for necessary information
    inputs = []
    for i in taglist:
        if i != None: 
            if i in ('squintmag', 'squintmag_uc', 'squintangle', 'squintangle_uc'):
                inputs.append('squintaz')
                inputs.append('squintel')
            elif i in ('feed'):
                inputs.append('round(feed,0) as feed')
            elif i in ('feedrev'):
                inputs.append('round(feed,1) as feedrev')
            elif i in ('antfeed'):
                inputs.append('1000*antnum + round(feed,0) as antfeed')
            elif i in ('antfeedrev'):
                inputs.append('1000*antnum + round(feed,1) as antfeedrev')
            else:
                inputs.append(i)
            
    if 'squintaz' in inputs: inputs.append('squintaz_uc')
    if 'squintel' in inputs: inputs.append('squintel_uc')
    if 'x_width_az' in inputs: inputs.append('x_width_az_uc')
    if 'x_width_el' in inputs: inputs.append('x_width_el_uc')
    if 'y_width_az' in inputs: inputs.append('y_width_az_uc')
    if 'y_width_el' in inputs: inputs.append('y_width_el_uc')
                
    # Get rid of duplicates
    inputs = list(set(inputs))
    
    # Take out flaggeg data
    if exclude_flagged:
        if wherecmd == '':
            wherecmd = 'WHERE flag=0'
        else:
            wherecmd = 'WHERE flag=0 AND (' + wherecmd[6:] + ')'
            
    sql_cmd = 'SELECT ' + ','.join(inputs) + ' FROM runs NATURAL JOIN obs ' + wherecmd
    cursor.execute(sql_cmd)
    
    # Turn into ndarray
    sqldata = sqlitecursor_to_ndarray(cursor)
    connection.close()
    
    # Create new array with derived tags, if needed
    getmag = 'squintmag' in taglist or 'squintmag_uc' in taglist
    getangle = 'squintangle' in taglist or 'squintangle_uc' in taglist

    extra_dtypes = []

    if getmag:
        extra_dtypes.append(('squintmag', '<f8'))
        extra_dtypes.append(('squintmag_uc', '<f8'))
    if getangle:
        extra_dtypes.append(('squintangle', '<f8'))
        extra_dtypes.append(('squintangle_uc', '<f8'))

    if len(extra_dtypes) == 0:
        return sqldata
    else:
        dtypes = eval(str(sqldata.dtype)) + extra_dtypes
        data = np.zeros(np.size(sqldata), dtype=dtypes)
        for i in sqldata.dtype.names:
            data[i] = sqldata[i]
        
    if getmag:
        data['squintmag'] = np.sqrt(data['squintaz'] ** 2 + data['squintel'] ** 2)
        data['squintmag_uc'] = np.sqrt((data['squintaz_uc'] * data['squintaz']) ** 2 + 
                                       (data['squintel_uc'] * data['squintel']) ** 2) / data['squintmag']

    if getangle:
        data['squintangle'] = np.arctan2(data['squintel'], data['squintaz'])
        el_over_az_uc = (np.sqrt((data['squintel_uc'] / data['squintel']) ** 2 +
                                 (data['squintaz_uc'] / data['squintaz']) ** 2)
                         * data['squintel'] / data['squintaz'])
        data['squintangle_uc'] = (el_over_az_uc /
                                  (1 + (data['squintel'] / data['squintaz']) ** 2))

    return data
    

def hexfromtxt(fname, dtype=float, names=None, skip_header=0, colnum=0):    
    """
    hexfromtxt
    ==========
    
    PURPOSE:
        Mimic numpy.genfromtxt from 1.4+ versions of numpy
        Read whitespace-separated text files into ndarray
        
    CALLING SEQUENCE:
        filearray = hexfromtxt(fname, dtype=float, names=None, skip_header=0, colnum=0)
        
    INPUTS:
        fname       :=  path to text file
        dtype       :=  list of column datatypes as strings, or string of datatypes separated by commas (no spaces)
                        (default all floats)
        names       :=  list of column names, or string of names separated by commas (no spaces)
                        (default 'f0','f1',...,'fn')
        skip_header :=  number of lines to skip at top of file
        colnum      :=  number of columns to look for
                        (if not specified, will look at dtype then names then first row after header)
    """

    # Track if any rows have wrong number of columns
    anybad = False

    # Read lines from file, skipping header
    file = open(fname, 'r')
    fileread = file.readlines()[skip_header:]
    file.close ()
    
    # Bail if empty
    if fileread == []: return None
    
    # Turn entries into lists of strings
    fileread = [i.strip().split() for i in fileread]
    
    # Change string dtype, names into list
    if type(dtype) == type(str()): dtype = dtype.split(',')
    if type(names) == type(str()): names = names.split(',')
    
    # If colnum not specified, take as number of dtype, then names, then first element
    if colnum == 0: 
        if type(dtype) == type(list()): colnum = len(dtype)
        elif type(names) == type(list()): colnum = len(names)
        else: colnum = len(fileread[0])
    
    # Get lists of proper-length elements
    filelist = []
    for i in fileread:
        if len(i) == colnum:
            filelist.append(i)
        else:
            anybad = True

    if anybad:
        print '   Some rows didn\'t have the right number of columns.'

    # Create proper ndarray
    if names == None:
        if dtype == float: passdtype = (',f' * colnum)[1:]
        else: passdtype = ','.join(dtype)
    elif dtype == float:
        passdtype = zip(names, (',f' * colnum)[1:].strip().split(','))
    else: passdtype = zip(names, dtype)
    
    filearray = np.zeros(len(filelist), dtype=passdtype)
    
    # Add data to array
    for i in xrange(len(filelist)):
        for j in xrange(colnum):
            filearray[i][j] = filelist[i][j]
    
    return filearray
