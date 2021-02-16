@Grapes([
    @Grab(group='org.semanticweb.elk', module='elk-owlapi', version='0.4.3'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='5.1.14'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='5.1.14'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='5.1.14'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-parsers', version='5.1.14'),
    @Grab(group='net.sourceforge.owlapi', module='owlapi-distribution', version='5.1.14'),
    @GrabConfig(systemClassLoader=true)
])

import org.semanticweb.owlapi.model.IRI
import org.semanticweb.owlapi.model.parameters.*
import org.semanticweb.elk.owlapi.*
import org.semanticweb.elk.reasoner.config.*
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.reasoner.*
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.io.*
import org.semanticweb.owlapi.owllink.*
import org.semanticweb.owlapi.util.*
import org.semanticweb.owlapi.search.*
import org.semanticweb.owlapi.manchestersyntax.renderer.*
import org.semanticweb.owlapi.reasoner.structural.*


def hpPath = '../../miesim/hp-inferred.owl'
def manager = OWLManager.createOWLOntologyManager()
def fac = manager.getOWLDataFactory()
def hpo = manager.loadOntologyFromOntologyDocument(new File(hpPath))

def allParents = { o, c ->
  o.getSubClassAxiomsForSubClass(fac.getOWLClass(IRI.create(c))).collect { a -> 
    def spc = a.getSuperClass()
    if(spc.getClassExpressionType().toString() =~ 'AllValuesFrom') {
      spc.getFiller().collect { fc -> fc.getIRI().toString() } 
    } else {
      spc.getIRI().toString() 
    }
  }.flatten()
}
def rewriteTermID = { l ->
  'http://purl.obolibrary.org/obo/' + l.replace(':', '_')
}

def class_parents = [:]
new File('/home/slater/miesim/dphens/sim_all_no_na.lst').splitEachLine(',') {
  def c1 = rewriteTermID(it[0])
  def c2 = rewriteTermID(it[1])
  def lcs

  if(!class_parents[c1]) { class_parents[c1] = [ [c1] ] }
  if(!class_parents[c2]) { class_parents[c2] = [ [c2] ] }

  // This will always terminate, since eventually the LCS will be owl:Thing
  def i = 1
  while(!lcs) {
    if(!class_parents[c1][i]) {
      class_parents[c1] << class_parents[c1][i-1].collect { c -> allParents(hpo, c) }.flatten()
    }
    if(!class_parents[c2][i]) {
      class_parents[c2] << class_parents[c2][i-1].collect { c -> allParents(hpo, c) }.flatten()
    }

    def all_c1 = (0..i).collect { r -> class_parents[c1][r] }.flatten()
    def all_c2 = (0..i).collect { r -> class_parents[c2][r] }.flatten()

    lcs = all_c1.find { ac -> all_c2.contains(ac) }
    i++
  }

  println "${it.join(',')},$lcs"
}

