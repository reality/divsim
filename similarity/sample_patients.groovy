// Here we will select 1,000 patient visits...

def all_patient_visits = [:]
def sampled_patients = []
def mapped = [:]
new File('../../miesim/doid_icd_mappings.txt').splitEachLine('\t') { mapped[it[0].replace('.','')] = it[1] }

new File('../../miesim/DIAGNOSES_ICD.csv').splitEachLine(',') { f ->
  def key = f[1] + '_' + f[2]
  if(!all_patient_visits.containsKey(key)) {
    all_patient_visits[key] = [] 
  }
  all_patient_visits[key] << f[4]
}

def rng = new Random()
pKeys = all_patient_visits.keySet().collect()

def already = []

while(sampled_patients.size() < 500) {
  def key = pKeys[rng.nextInt(all_patient_visits.size())]
  def pt = all_patient_visits[key]
  if(!pt[0]) { continue; }
  def pDiag = pt[0].replaceAll('"', '')
  if(new File('../../miesim/new_texts/' + key + '.txt').exists() && mapped[pDiag] && !already.contains(key)) {
    sampled_patients << key + '\t' + pt[0] + '\t' + mapped[pDiag]
    already << key
  }
}

new File('sampled_patient_visits.csv').text = sampled_patients.join('\n')
