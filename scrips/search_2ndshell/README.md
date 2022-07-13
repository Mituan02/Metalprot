# This is the first method to search 2nd shell for metal binding, which is currently different from the current method using COMBS database.

# This method is more metal specific as the database is extracted from metal only database.

# To use this method, please check ../database_gneerate_2ndshell.

# It is important to analyze: 
  Is the metal binding need 2ndshell?  Yes, every metal binding has at least one 2nd shell contact.
  Is the metal binding 2ndshell from neighbor aa or remote aa? About 1/3 are within 7 aas.
  Is the metal specific 2nshell different from the general hbond? (Need to compare the two different database)
  Is the current COMBS database good enough for searching the metal binding 2ndshell? (The COMBS database remove neighbor aa contact)