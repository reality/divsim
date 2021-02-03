def keepFile = args[0]
def keep = []
new File(keepFile).splitEachLine('\t') {
  keep << it[0]
}
new File('./dphens_merged.txt').splitEachLine('\t') { profile ->
  println profile[0] + '\t' + profile[1].split(',').findAll { keep.contains(it) }.join(',')
}
