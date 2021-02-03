def profiles = [:]
new File('./dphens_baseline.txt').splitEachLine('\t') {
  profiles[it[0]] = it[1].tokenize(',')
}

new File('./similarity/trainpat/patient_disease_profiles.lst').splitEachLine('\t') {
  if(!profiles.containsKey(it[0])) { profiles[it[0]] = [] }
  profiles[it[0]] += it[1].tokenize(',')
  profiles[it[0]].unique(true)
}

new File('./dphens_with_patient_training.txt').text = profiles.collect { k, v ->
  "$k\t${v.join(',')}"
}.join('\n')
