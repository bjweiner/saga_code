#!python
#
# Take a list of SAGA reduced data files and extract information like
# objid and fiber number so we can find spectra, cross-ref spectra
# to catalog, etc.  The output will be an astropy Table written as
# a FITS table.
#
# The fields extracted for each target are name, ra, dec, mag,
#  filename, aperture (ie number within the file)
#
# BJW, Jan-Mar 2018

from astropy.io import fits
# from astropy.io import ascii
from astropy.table import Table, vstack
import os

# read a list of filenames from a file with one name per line
def read_files(fname):
    listfile = open(fname,'r')
    flist = listfile.read().splitlines()
    listfile.close()
    return flist

# read two list of filenames from a file with two names per line
def read_two_files(fname):
    listfile = open(fname,'r')
    flist = []
    clist = []
    for line in listfile:
        fields = (line.strip()).split()
        if len(fields) < 2:
            print 'I didnt get two filenames from ', fname
        flist.append(fields[0])
        clist.append(fields[1])
    listfile.close()
    return flist, clist

# AAT files like OBrother_1_MG.fits.gz have
# ext 1 = 4940 x 400 image, spectra x object
# ext 2 = 17 x 400 table, catdata
# get fields name (string), ra, dec, magnitude
def getstruct_spectra_metadata_aat(fname):
    hdulist = fits.open(fname)
    header = hdulist[1].header
    cat = hdulist[2]
    outstruct = []
    ncat = len(cat.data)
    for i in range(ncat):
        cat1 = cat.data[i]
        struct1 = {'name':cat1['name'], 'ra':cat1['ra'], 'dec':cat1['dec'], 'mag':cat1['magnitude'], 'filename':fname, 'aperture':i}
        #struct1 = {'name':'', 'ra':0.0, 'dec':0.0, 'mag':0.0, 'filename':'', 'aperture':0}
        outstruct.append(struct1)
    hdulist.close()
    return outstruct

# Create an astropy table suitable for writing into a fits table
def get_spectra_metadata_aat(fname):
    hdulist = fits.open(fname)
    header = hdulist[1].header
    cat = hdulist[2]
    ncat = len(cat.data)
    # Make a list ncat long with filename in each entry
    fname_list = [fname] * ncat
    aper_list = range(ncat)
    table_aat = Table([cat.data['name'], cat.data['ra'], cat.data['dec'], cat.data['magnitude'], fname_list, aper_list], names=('name', 'ra', 'dec', 'mag', 'filename', 'aperture'))
    hdulist.close()
    print table_aat[0]
    return table_aat

# MMT files like hecto_2017a/2017.0224/reduction/0000/spHect-2017.0224_1.fits
# ext 1 = 4608 x 300 image, spectra x object
# ext 2,3,4 = also images
# ext 5 = 18 x 300 table, catdata
# but unfortunately the table does not include the object ID!
# get fields name (string), ra, dec, magnitude
# need to get names from another catalog. The spZall*.fits doesn't have
# the names and has many more entries, like 51 x 300.
# The plugdir/ScoobyDoo_jan2017_1_cat  does have the object names
def getstruct_spectra_metadata_mmt_noname(fname):
    hdulist = fits.open(fname)
    header = hdulist[1].header
    cat = hdulist[5]
    outstruct = []
    ncat = len(cat.data)
    for i in range(ncat):
        cat1 = cat.data[i]
        struct1 = {'name':'', 'ra':cat1['ra'], 'dec':cat1['dec'], 'mag':cat1['rmag'], 'filename':fname, 'aperture':i}
        #struct1 = {'name':'', 'ra':0.0, 'dec':0.0, 'mag':0.0, 'filename':'', 'aperture':0}
        outstruct.append(struct1)
    hdulist.close()
    return outstruct

# The plugdir/ScoobyDoo_jan2017_1_cat  does have the object names
# Here we also read the catalog text file to get those names.
# Its columns are: aperture, ra-hours, dec, mag, fibernum?, iftarget, objname
def getstruct_spectra_metadata_mmt_cat(fname,catname=''):
    hdulist = fits.open(fname)
    header = hdulist[1].header
    cat = hdulist[5]
    outstruct = []
    ncat = len(cat.data)
    if catname != '':
        f = open(catname, 'r')
        objnames = []
        mags = []
        for line in f:
            fields = (line.strip()).split()
            mag1 = float(fields[3])
            # Trap lines with a blank name (unused fibers)
            if len(fields) >= 7:
                name1 = fields[6]
            else:
                name1 = 'unused'
            objnames.append(name1)
            mags.append(mag1)
        f.close()
    else:
        print 'Need an MMT catalog name to get objIDs for ',fname
    if ncat != len(objnames):
        print 'Lengths of spectrum file and catalog do not match', fname
        
    for i in range(ncat):
        cat1 = cat.data[i]
        struct1 = {'name':objnames[i], 'ra':cat1['ra'], 'dec':cat1['dec'], 'mag':cat1['rmag'], 'filename':fname, 'aperture':i}
        #struct1 = {'name':'', 'ra':0.0, 'dec':0.0, 'mag':0.0, 'filename':'', 'aperture':0}
        outstruct.append(struct1)
    hdulist.close()
    return outstruct

# create a table for MMT, need to read two files
# read the catalog text file to get those names.
# Its columns are: aperture, ra-hours, dec, mag, fibernum?, iftarget, objname

def get_spectra_metadata_mmt_cat(fname,catname=''):
    hdulist = fits.open(fname)
    header = hdulist[1].header
    cat = hdulist[5]
    ncat = len(cat.data)
    # Make a list ncat long with filename in each entry
    fname_list = [fname] * ncat
    aper_list = range(ncat)
    if catname != '':
        f = open(catname, 'r')
        objnames = []
        mags = []
        for line in f:
            fields = (line.strip()).split()
            mag1 = float(fields[3])
            # Trap lines with a blank name (unused fibers)
            if len(fields) >= 7:
                name1 = fields[6]
            else:
                name1 = 'unused'
            objnames.append(name1)
            mags.append(mag1)
        f.close()
    else:
        print 'Need an MMT catalog name to get objIDs for ',fname
    if ncat != len(objnames):
        print 'Lengths of spectrum file and catalog do not match', fname

    table_mmt = Table([objnames, cat.data['ra'], cat.data['dec'], mags, fname_list, aper_list], names=('name', 'ra', 'dec', 'mag', 'filename', 'aperture'))
    hdulist.close()
    print table_mmt[0]
    return table_mmt
  


####

# write the data structure / table to a file as a fits table

def write_spectra_cat(struct,outname):
    # if struct is an astropy table
    struct.write(outname, format='fits')
    return

def do_aat_files(lname):
    afilelist = read_files(lname)
    # this made a list, but I want a fits table
    # outtable = []
    ifirst = 1
    for fname in afilelist:
        table1 = get_spectra_metadata_aat(fname)
        print "Read AAT rows ",len(table1)
        if ifirst == 1:
            outtable = table1
            ifirst = 0
        else:
            outtable = vstack([outtable,table1])
    print outtable[0]        
    return outtable

def do_mmt_files(lname):
    mfilelist, mcatlist = read_two_files(lname)
    # this made a list, but I want a fits table
    # outtable = []
    ifirst = 1
    for fname, cname in zip(mfilelist, mcatlist):
        table1 = get_spectra_metadata_mmt_cat(fname,catname=cname)
        print "Read MMT rows ",len(table1)
        if ifirst == 1:
            outtable = table1
            ifirst = 0
        else:
            outtable = vstack([outtable,table1])
    print outtable[0]
    return outtable

def main():
    fname1 = raw_input('Enter file with list of AAT files: ')
    fname2 = raw_input('Enter file with list of MMT spec and cat files: ')
    outname = raw_input('Name for output fits file: ')
    aattable = do_aat_files(fname1)
    mmttable = do_mmt_files(fname2)
    print "AAT length ",len(aattable),"  MMT length ",len(mmttable)
    outtable = vstack([aattable,mmttable])
    print "total length ",len(outtable)
    print outtable[0]
    print outtable[len(outtable)-1]
    if outname =='':
        outname = 'spectrum_table.fits'
    try:
        outtable.write(outname, format='fits')
    except:
        outtable.write("tmp_"+outname, format='fits')
        
    
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
