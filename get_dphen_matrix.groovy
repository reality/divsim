@Grab(group='org.semanticweb.elk', module='elk-owlapi', version='0.4.3')
@Grab(group='net.sourceforge.owlapi', module='owlapi-api', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-apibinding', version='4.2.5')
@Grab(group='net.sourceforge.owlapi', module='owlapi-impl', version='4.2.5')
@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')

import groovyx.gpars.*
import org.codehaus.gpars.*
import java.util.concurrent.*
 
import slib.graph.algo.utils.GAction;
import slib.graph.algo.utils.GActionType;

import org.openrdf.model.URI;
import slib.graph.algo.accessor.GraphAccessor;
//import slib.graph.algo.extraction.utils.GAction;
//import slib.graph.algo.extraction.utils.GActionType;
import slib.graph.algo.validator.dag.ValidatorDAG;
import slib.graph.io.conf.GDataConf;
import slib.graph.io.conf.GraphConf;
import slib.graph.io.loader.GraphLoaderGeneric;
import slib.graph.io.util.GFormat;
import slib.graph.model.graph.G;
import slib.graph.model.impl.graph.memory.GraphMemory;
import slib.graph.model.impl.repo.URIFactoryMemory;
import slib.graph.model.repo.URIFactory;
import slib.sml.sm.core.engine.SM_Engine;
import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Topo;
import slib.sml.sm.core.metrics.ic.utils.ICconf;
import slib.graph.algo.extraction.utils.*
import slib.sglib.io.loader.*
import slib.sml.sm.core.metrics.ic.utils.*
import slib.sml.sm.core.utils.*
import slib.sglib.io.loader.bio.obo.*
import org.openrdf.model.URI
import slib.graph.algo.extraction.rvf.instances.*
import slib.sglib.algo.graph.utils.*
import slib.utils.impl.Timer
import slib.graph.algo.extraction.utils.*
import slib.graph.model.graph.*
import slib.graph.model.repo.*
import slib.graph.model.impl.graph.memory.*
import slib.sml.sm.core.engine.*
import slib.graph.io.conf.*
import slib.graph.model.impl.graph.elements.*
import slib.graph.algo.extraction.rvf.instances.impl.*
import slib.graph.model.impl.repo.*
import slib.graph.io.util.*
import slib.graph.io.loader.*

import slib.sml.sm.core.metrics.ic.utils.IC_Conf_Corpus;
import slib.sml.sm.core.utils.SMConstants;
import slib.sml.sm.core.utils.SMconf;
import slib.utils.ex.SLIB_Exception;
import slib.utils.impl.Timer;

def oFile = args[0]

def patient_visit_diagnoses = [:]
def pcounts = [:]

ConcurrentHashMap aMap = [:]

def cList = []

def dMap = [:]
new File('./dphens_merged.txt').splitEachLine('\t') {
  dMap[it[0]] = [] 
  it[1].split(',').each { p ->
    p = p.replace('<', '').replace('>', '').tokenize('/').last()
    def z = p.tokenize('_')
    p = z[0]+':' + z[1]

    dMap[it[0]] << p
    if(!cList.contains(p)) {
      cList << p 
    }
  }
}

println 'writing the annotation file now'
def sWriter = new BufferedWriter(new FileWriter('annot.tsv'))
def oo = ''
def y = 0

dMap.each { a, b ->
  println "(${++y}/${aMap.size()})"
  if(!b.any{ it.indexOf('txt') != -1}) {
    sWriter.write('http://reality.rehab/ptvis/' + a + '\t' + b.join(';') + '\n')
  }
}
sWriter.flush()
sWriter.close()
println 'done'

URIFactory factory = URIFactoryMemory.getSingleton()

def ontoFile = oFile
def graphURI = factory.getURI('http://HP/')
factory.loadNamespacePrefix("HP", graphURI.toString());

G graph = new GraphMemory(graphURI)

def dataConf = new GDataConf(GFormat.RDF_XML, ontoFile)
def actionRerootConf = new GAction(GActionType.REROOTING)
actionRerootConf.addParameter("root_uri", "http://purl.obolibrary.org/obo/HP_0000001"); // phenotypic abnormality
//actionRerootConf.addParameter("root_uri", "DOID:4"); // phenotypic abnormality

def gConf = new GraphConf()
gConf.addGDataConf(dataConf)
gConf.addGAction(actionRerootConf)
//def gConf = new GraphConf()
//gConf.addGDataConf(dataConf)
def annot = 'annot.tsv'
gConf.addGDataConf(new GDataConf(GFormat.TSV_ANNOT, annot));

GraphLoaderGeneric.load(gConf, graph)

println graph.toString()

def roots = new ValidatorDAG().getTaxonomicRoots(graph)
println roots

//def icConf = new IC_Conf_Corpus(SMConstants.FLAG_IC_ANNOT_RESNIK_1995)
def icConf = new IC_Conf_Topo(SMConstants.FLAG_ICI_ZHOU_2008)
//def icConf = new IC_Conf_Topo(SMConstants.FLAG_ICI_SANCHEZ_2011)
def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_RESNIK_1995, icConf)
//def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_LIN_1998, icConf)
//def smConfGroupwise = new SMconf(SMConstants.FLAG_SIM_GROUPWISE_AVERAGE, icConf)
def smConfGroupwise = new SMconf(SMConstants.FLAG_SIM_GROUPWISE_BMA, icConf)
// FLAG_SIM_GROUPWISE_AVERAGE_NORMALIZED_GOSIM

//def smConfPairwise = new SMconf(SMConstants.FLAG_SIM_PAIRWISE_DAG_NODE_JIANG_CONRATH_1997_NORM , icConf)


def out = []
def z = 0

def outWriter = new BufferedWriter(new FileWriter('sim.lst'), 1024 * 1024);

def engine = new SM_Engine(graph)

def getURIfromTerm = { term ->
  term = term.tokenize(':')
  return factory.getURI("http://purl.obolibrary.org/obo/HP_" + term[1])
}


cList.eachWithIndex { a, ii ->
println "${++z}/${cList.size()}"
  def aList = []
  cList.eachWithIndex { b, j ->
    if(ii > j) { return; }
    def value = 'NA'
    try {
    value = engine.compare(smConfPairwise, getURIfromTerm(a), getURIfromTerm(b))
    } catch(e) { println 'err comp'; value = 'NA'; }
    aList << [a, b, value] 
  }

  aList = aList.toSorted { 
    def c = 0
    if(it[2] != 'NA') { c = it[2] } 
    c 
  }.reverse() //[0..10]

  aList.eachWithIndex { it, i -> 
    outWriter.write(it.join(',') + ',' + (i+1) + '\n')
  }
}

outWriter.flush()
outWriter.close()

println 'done'
