import groovyx.gpars.*
import org.codehaus.gpars.*
import java.util.concurrent.*

def sim = [:]
new File('pre_weights.tsv').splitEachLine('\t') {
  if(it[6] == 'similarity') { return; } // yes this is dumb, it's 1am
  if(!sim.containsKey(it[0])) {
    sim[it[0]] = [:] 
  }
  def value = it[6].toFloat()
  sim[it[0]][it[1]] = value
  if(!sim.containsKey(it[1])) {
    sim[it[1]] = [:] 
  }
  sim[it[1]][it[0]] = value
}

println 'Loaded weights.'

def patients = [:]
new File('./annotations.txt').splitEachLine('\t') {
  def uid = it[0].tokenize('.')[0]
  if(sim[it[1]]) {
    if(!patients.containsKey(uid)) {
      patients[uid] = []
    }
    patients[uid] << it[1]
  }
}

println 'Loaded patients.'

def diagnoses = [:]
new File('sampled_patient_visits.csv').splitEachLine('\t') {
  diagnoses[it[0]] = it[2]
}

println 'Loaded patient diagnoses.'

def diseases = [:]
new File('../dphens_with_patient_training.txt').splitEachLine('\t') {
  if(it[1]) { 
  diseases[it[0]] = it[1].tokenize(',')
  }
}

println 'Loaded diseases.'

def outWriter = new BufferedWriter(new FileWriter('new_test_matrix_with_training.lst'), 1024 * 1024)

def z = 0
GParsPool.withPool(45) { p ->  // im a monster
patients.eachParallel { pid, u1 ->
  println "(${++z}/${patients.size()})"

  def aList = []
  diseases.each { doid, u2 ->
    aList << [
      doid,
      ((u1.inject(0) { sum, uri ->
        sum += u2.collect { uri2 -> 
          if(sim[uri] && sim[uri][uri2]) {
            sim[uri][uri2]
          } else { 0} }.max()
      } + u2.inject(0) { sum, uri ->
        sum += u1.collect { uri2 -> 
          if(sim[uri] && sim[uri][uri2]) {
            sim[uri][uri2]
          } else { 0 }}.max()
      }) / (u1.size() + u2.size())) // /
      // there's some extra effect with NAs here; bad
    ]
  }

  aList.toSorted { it[1] }.reverse().eachWithIndex { it, i ->
    def match = false
    if(diagnoses[pid] == it[0]) {
      match = true 
    }

    outWriter.write(pid + ',' + it[0] + ',' + (i+1) + ','+ it[1]+ ',' + match + '\n')
  }
}
}

outWriter.flush()
outWriter.close()

println 'Done.'
