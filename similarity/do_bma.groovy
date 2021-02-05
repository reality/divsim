def sim = [:]
new File('pre_weights.tsv').splitEachLine('\t') {
  if(!sim.containsKey(it[0])) {
    sim[it[0]] = [:] 
  }
  sim[it[0]][it[1]] = it[6] 
  if(!sim.containsKey(it[1])) {
    sim[it[1]] = [:] 
  }
  sim[it[1]][it[0]] = it[6] 
}

println 'Loaded weights.'

def patients = [:]
new File('annotations.tsv').splitEachLine('\t') {
  def uid = it[0].tokenize('.')[0]
  if(!patients.containsKey(uid)) {
    patients[uid] = []
  }
  patients[uid] << it[1]
}

println 'Loaded patients.'

def diagnoses = [:]
new File('sampled_patient_visits.csv').splitEachLine('\t') {
  diagnoses[it[0]] = it[2]
}

println 'Loaded patient diagnoses.'

def diseases = [:]
new File('../dphens_baseline.txt').splitEachLine('\t') {
  diseases[it[0]] = it[1].tokenize(',')
}

println 'Loaded diseases.'

GParsPool.withPool(45) { p ->  // im a monster
patients.eachParallel { pid, u1 ->
  println "(${++z}/${aMap.size()})"

  def aList = []
  diseases.each { doid, u2 ->
    aList << [
      doid,
      ((u1.inject(0) { sum, uri ->
        sum += u2.collect { uri2 -> cSim[uri][uri2] }.max()
      } + u2.inject(0) { sum, uri ->
        sum += u1.collect { uri2 -> cSim[uri][uri2] }.max()
      }) / (u1.size() + u2.size())) // /
    ]
  }

  def match = false
  if(diagnoses[pid] == doid) {
    match = true 
  }

  aList.toSorted { it[1] }.reverse().eachWithIndex { i, it ->
    outWriter.write(g1 + ',' + it[0] + ',' + (i+1) + ','+ it[1]+ ',' + match + '\n')
  }
}
}

println 'Done.'
