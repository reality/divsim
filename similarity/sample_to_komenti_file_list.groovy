def out = []
new File('./sampled_patient_visits.csv').splitEachLine('\t') {
  out << '../../miesim/new_texts/' + it[0] + '.txt'
}

new File('annotate_flist.txt').text = out.join('\n')
