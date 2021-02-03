def chars = [:]
def trueBeforeTen = [:]
def scores = []
new File(args[0]).splitEachLine(','){
  if(!chars.containsKey(it[0])) {
    chars[it[0]] = 0
    trueBeforeTen[it[0]] = 0
  }
  if(it[3] == 'true' && chars[it[0]] == 0) {
    chars[it[0]] = 1/it[4].toInteger()
  } //
  if(it[4].toInteger() <= 10) {
    if(it[3] == 'true')  {
      trueBeforeTen[it[0]] = 1
    }
  }
  scores << it[2].toFloat()
}
println scores.sum() / scores.size() //

println chars.collect {k,v->v}.sum() / chars.size() //
def o = trueBeforeTen.collect {k,v->v}.sum() / chars.size() //
println o*100
