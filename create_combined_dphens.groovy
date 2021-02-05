def profiles = [:]
new File('./dphens_merged_no7.5p_belowmean.lst').splitEachLine('\t') {
  if(!it[1]) {
 profiles[it[0]]  = []
  } else{
  profiles[it[0]] = it[1].tokenize(',')
  }
}

new File(args[0]).splitEachLine('\t') {
  if(!profiles.containsKey(it[0])) { profiles[it[0]] = [] }
  profiles[it[0]] += it[1].tokenize(',')
  profiles[it[0]].unique(true)
}

new File(args[1]).text = profiles.collect { k, v ->
  "$k\t${v.join(',')}"
}.join('\n')
