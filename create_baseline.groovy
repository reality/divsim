def keepFile = './nona_list.txt'
def keep = []
new File(keepFile).splitEachLine(',') {
  keep << 'http://purl.obolibrary.org/obo/'+it[0].replace(':','_')
}
new File('./dphens_merged.txt').splitEachLine('\t') { profile ->
  println profile[0] + '\t' + profile[1].split(',').findAll { keep.contains(it) }.join(',')
}
