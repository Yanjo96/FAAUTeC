from pylogeny.executable import executable, aTemporaryDirectory, treepuzzle, E_TREEPUZZ
from pylogeny import tree
from os.path import abspath, isdir, isfile

class raxml_2(executable):

    ''' Denotes a single run of the RAxML executable. See
    http://sco.h-its.org/exelixis/software.html for more information on RAxML.
    Requires RAxML to be installed.
    For our purpose we need some extra functions which are not implemented in
    the original pylogeny code
    Therefore we added a raxml_2 class which can handle the constraint tree'''

    exeName = 'raxmlHPC'

    def __init__(self,inp_align,out_file,model=None,
                 is_Protein=True,interTrees=False,
                 alg=None,startingTree=None, grouping_constraint=None,
                 bipartition_filename=None, parsimony_seed=None,rapid=False,
                 slow=False, optimizeBootstrap=False,numboot=100, log=None,
                 wdir="."):
        self.alignment            = inp_align
        self.out                  = out_file
        self.model                = model
        self.protein              = is_Protein
        self.intermediates        = interTrees
        self.alg                  = alg
        self.startingTree         = startingTree
        self.grouping_constraint  = grouping_constraint
        self.bipartition_filename = bipartition_filename
        self.parsimony_seed       = parsimony_seed
        self.slowBoot             = slow
        self.rapidBoot            = rapid
        self.numBoot              = numboot
        self.optBoot              = optimizeBootstrap
        self.log                  = log
        self.workdir              = abspath(wdir)
        if not self.model:
            if is_Protein: self.model = 'PROTGAMMAJTT'
            else:          self.model = 'GTRGAMMA'

    def getInstructionString(self):
        s = '%s -s %s -n %s -m %s' %\
            (self.exeName,self.alignment,self.out,self.model)
        if self.alg: s += ' -f %s' % (self.alg)
        if self.intermediates: s += ' -j'
        if self.startingTree: s += ' -t %s' % (self.startingTree)
        if self.grouping_constraint: s += ' -g %s' % (self.grouping_constraint)
        if self.bipartition_filename: s += ' -z %s' % (self.bipartition_filename)
        if self.parsimony_seed: s += ' -p %s' % (self.parsimony_seed)
        if self.rapidBoot or self.slowBoot:
            if self.rapidBoot:  s += ' -x 234534251 -N %d' % (self.numBoot)
            elif self.slowBoot: s += ' -b 234534251 -N %d' % (self.numBoot)
        if self.optBoot: s += ' -k'
        if self.workdir: s += ' -w %s' % (self.workdir)
        if self.log: s += ' > %s' % (self.log)
        return s

    def runFunction(self):
        self.run()
        return self.getInstructionString()

class consel_2(executable):

    ''' Denotes a single run of the CONSEL workflow in order to acquire a
    confidence interval and perform an AU test on a set of trees. Requires
    CONSEL to be installed. '''

    def __init__(self,treeset,alignment, sitelh, name, consel_path):

        self.treeset     = treeset
        self.name        = name
        self.consel_path = consel_path
        self.alignment   = alignment
        self._out        = None
        self.sitelh      = sitelh
        self.raw         = None
        self.rmt         = None
        self.pv          = None
        self.auvals      = list()
        self.interval    = tree.treeSet()
        self.rejected    = tree.treeSet()
        self.instruction = ''

    def getInstructionString(self):

        ''' Get the instruction string.

        :return: a string (of a UNIX command)

        '''

        return self.instruction

    def _convertRawData(self):

        if not isfile('%s.mt' % (self.name)):
            self.instruction = '%s/seqmt --puzzle %s %s.mt' %\
                (self.consel_path,self.sitelh,self.name)
            self._out = self.run()
            if not isfile('%s.mt' % (self.name)):
                raise IOError("CONSEL 'seqmt' did not create mt file.")
        self.raw = '%s.mt' % (self.name)

    def _createReplicates(self):
        if not isfile('%s.rmt' % (self.name)):
            self.instruction = '%s/makermt %s' % (self.consel_path, self.raw)
            self._out = self.run()
            if not isfile('%s.rmt' % (self.name)):
                raise IOError("CONSEL 'markermt' did not create rmt file.")
        self.rmt = '%s.rmt' % (self.name)

    def _run(self):

        if not isfile('%s.pv' % (self.name)):
            self.instruction = '%s/consel %s.rmt' % (self.consel_path,self.name)
            self._out = self.run()
        self.pv = '%s.pv' % (self.name)

    def _getAU(self):

        self.instruction = '%s/catpv %s' % (self.consel_path, self.pv)
        pvout = self.run()
        for line in pvout.split('\n'):
            spl = line.split()
            if len(spl) < 3: continue
            elif not spl[2].isdigit(): continue
            it = (int(spl[2])-1,float(spl[4]))
            self.auvals.append(it)
            if it[1] >= 0.05: self.interval.addTree(self.treeset[it[0]])
            else: self.rejected.addTree(self.treeset[it[0]])
        with open(self.name + ".consel", "w") as consel_out:
            for line in pvout:
                consel_out.write(line)
        self._out = pvout

    def getRejected(self):

        ''' If an AU test has already been performed, return the set of trees
        that were rejected by the test.

        :return: a :class:`.tree.treeSet` object or None if no test was done yet

        '''

        if (self._out != None):
            return self.rejected
        return None

    def getInterval(self):

        ''' Compute the AU test. Return the interval of trees as a tree set.

        :return: a :class:`.tree.treeSet` object

        '''
        if not isfile('%s.trees.sitelh' % (self.name)):
            self.treefile = self.treeset.toTreeFile(self.name + '.trees')
            self.sitelh = treepuzzle(self.alignment,self.treefile).getSiteLikelihoodFile()
        else: self.sitelh = '%s.trees.sitelh' % (self.name)
        self._convertRawData()
        self._createReplicates()
        self._run()
        self._getAU()
        return self.interval

    def getLog(self):
        log = []
        self._convertRawData()
        log.append(self.instruction)
        self._createReplicates()
        log.append(self.instruction)
        self._run()
        log.append(self.instruction)
        self._getAU()
        log.append(self.instruction)
        return log
