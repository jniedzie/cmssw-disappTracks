
COM_ENERGY = 13000.
CROSS_SECTION = 1.0
MCHI = 300  # GeV
CTAU = 100  # mm
SLHA_TABLE="""
  #  ISAJET SUSY parameters in SUSY Les Houches Accord 2 format
  #  Created by ISALHA 2.0 Last revision: C. Balazs 21 Apr 2009
  Block SPINFO   # Program information
  1   ISASUGRA from ISAJET          # Spectrum Calculator
  2   7.80   29-OCT-2009 12:50:36   # Version number
  Block MODSEL   # Model selection
  1     3   # Minimal anomaly mediated (AMSB) model
  Block SMINPUTS   # Standard Model inputs
  1     1.27843719E+02   # alpha_em^(-1)
  2     1.16570000E-05   # G_Fermi
  3     1.17200002E-01   # alpha_s(M_Z)
  4     9.11699982E+01   # m_{Z}(pole)
  5     4.19999981E+00   # m_{b}(m_{b})
  6     1.73070007E+02   # m_{top}(pole)
  7     1.77699995E+00   # m_{tau}(pole)
  Block MINPAR   # SUSY breaking input parameters
  1     1.50000000E+03   # m_0
  2     1.03650000E+05   # m_{3/2}
  3     5.00000000E+00   # tan(beta)
  4     1.00000000E+00   # sign(mu)
  Block EXTPAR   # Non-universal SUSY breaking parameters
  0     1.21940317E+16   # Input scale
  Block MASS   # Scalar and gaugino mass spectrum
  #  PDG code   mass                 particle
  24     8.04229965E+01   #  W^+
  25     1.13281731E+02   #  h^0
  35     2.31703979E+03   #  H^0
  36     2.30156934E+03   #  A^0
  37     2.31312744E+03   #  H^+
  1000001     2.49383179E+03   #  dnl
  1000002     2.49261523E+03   #  upl
  1000003     2.49383179E+03   #  stl
  1000004     2.49261548E+03   #  chl
  1000005     2.15123779E+03   #  b1
  1000006     1.73482275E+03   #  t1
  1000011     1.43890710E+03   #  el-
  1000012     1.43549438E+03   #  nuel
  1000013     1.43890710E+03   #  mul-
  1000014     1.43549438E+03   #  numl
  1000015     1.42184009E+03   #  tau1
  1000016     1.43293103E+03   #  nutl
  1000021     2.18962207E+03   #  glss
  1000022     2.99864838E+02   #  z1ss
  1000023     9.52063293E+02   #  z2ss
  1000024     3.00137323E+02   #  w1ss
  1000025    -1.74983081E+03   #  z3ss
  1000035     1.75280737E+03   #  z4ss
  1000037     1.75418726E+03   #  w2ss
  2000001     2.52172412E+03   #  dnr
  2000002     2.50711743E+03   #  upr
  2000003     2.52172412E+03   #  str
  2000004     2.50711792E+03   #  chr
  2000005     2.50527832E+03   #  b2
  2000006     2.17448901E+03   #  t2
  2000011     1.42790625E+03   #  er-
  2000013     1.42790625E+03   #  mur-
  2000015     1.43795154E+03   #  tau2
  Block ALPHA   # Effective Higgs mixing parameter
  -1.98217750E-01   # alpha
  Block STOPMIX   # stop mixing matrix
  1  1     1.15278468E-01   # O_{11}
  1  2    -9.93333220E-01   # O_{12}
  2  1     9.93333220E-01   # O_{21}
  2  2     1.15278468E-01   # O_{22}
  Block SBOTMIX   # sbottom mixing matrix
  1  1     9.99966919E-01   # O_{11}
  1  2     8.13260302E-03   # O_{12}
  2  1    -8.13260302E-03   # O_{21}
  2  2     9.99966919E-01   # O_{22}
  Block STAUMIX   # stau mixing matrix
  1  1     3.27634096E-01   # O_{11}
  1  2     9.44804668E-01   # O_{12}
  2  1    -9.44804668E-01   # O_{21}
  2  2     3.27634096E-01   # O_{22}
  Block NMIX   # neutralino mixing matrix
  1  1    -1.78793352E-03   #
  1  2     9.98812914E-01   #
  1  3    -4.60493490E-02   #
  1  4     1.57782547E-02   #
  2  1     9.99014318E-01   #
  2  2     3.87472054E-03   #
  2  3     3.69140282E-02   #
  2  4    -2.43423060E-02   #
  3  1     8.95513594E-03   #
  3  2    -2.13845931E-02   #
  3  3    -7.06530809E-01   #
  3  4    -7.07302272E-01   #
  4  1     4.34374809E-02   #
  4  2    -4.35934253E-02   #
  4  3    -7.05217123E-01   #
  4  4     7.06315577E-01   #
  Block UMIX   # chargino U mixing matrix
  1  1    -9.97862458E-01   # U_{11}
  1  2     6.53490350E-02   # U_{12}
  2  1    -6.53490350E-02   # U_{21}
  2  2    -9.97862458E-01   # U_{22}
  Block VMIX   # chargino V mixing matrix
  1  1    -9.99705255E-01   # V_{11}
  1  2     2.42768954E-02   # V_{12}
  2  1    -2.42768954E-02   # V_{21}
  2  2    -9.99705255E-01   # V_{22}
  Block GAUGE Q=  1.85234131E+03   #
  1     3.57492119E-01   # g`
  2     6.52495980E-01   # g_2
  3     1.22099257E+00   # g_3
  Block YU Q=  1.85234131E+03   #
  3  3     8.58246148E-01   # y_t
  Block YD Q=  1.85234131E+03   #
  3  3     6.76042885E-02   # y_b
  Block YE Q=  1.85234131E+03   #
  3  3     5.14520593E-02   # y_tau
  Block HMIX Q=  1.85234131E+03   # Higgs mixing parameters
  1     1.74554053E+03   # mu(Q)
  2     5.00000000E+00   # tan(beta)(M_GUT)
  3     2.51296097E+02   # Higgs vev at Q
  4     5.29722150E+06   # m_A^2(Q)
  Block MSOFT Q=  1.85234131E+03   # DRbar SUSY breaking parameters
  1     9.59482605E+02   # M_1(Q)
  2     2.83875336E+02   # M_2(Q)
  3    -2.03255798E+03   # M_3(Q)
  31     1.43290662E+03   # MeL(Q)
  32     1.43290662E+03   # MmuL(Q)
  33     1.43040063E+03   # MtauL(Q)
  34     1.42463135E+03   # MeR(Q)
  35     1.42463135E+03   # MmuR(Q)
  36     1.41936206E+03   # MtauR(Q)
  41     2.38789746E+03   # MqL1(Q)
  42     2.38789746E+03   # MqL2(Q)
  43     2.05947778E+03   # MqL3(Q)
  44     2.40229883E+03   # MuR(Q)
  45     2.40229883E+03   # McR(Q)
  46     1.66603821E+03   # MtR(Q)
  47     2.41674878E+03   # MdR(Q)
  48     2.41674878E+03   # MsR(Q)
  49     2.42613452E+03   # MbR(Q)
  Block AU Q=  1.85234131E+03   #
  1  1     1.71933716E+03   # A_u
  2  2     1.71933716E+03   # A_c
  3  3     1.71933716E+03   # A_t
  Block AD Q=  1.85234131E+03   #
  1  1     4.09207031E+03   # A_d
  2  2     4.09207031E+03   # A_s
  3  3     4.09207031E+03   # A_b
  Block AE Q=  1.85234131E+03   #
  1  1     1.09036902E+03   # A_e
  2  2     1.09036902E+03   # A_mu
  3  3     1.09036902E+03   # A_tau
  #
  #
  #
  #                             =================
  #                             |The decay table|
  #                             =================
  #
  #         PDG            Width
  DECAY   1000024     %.9g # chargino decay
  #
  #
  """ % (1.97326979e-13 / CTAU)

import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *


generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(-1),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         SLHATableForPythia8 = cms.string('%s' % SLHA_TABLE),
                         comEnergy = cms.double(COM_ENERGY),
                         crossSection = cms.untracked.double(CROSS_SECTION),
                         maxEventsToPrint = cms.untracked.int32(0),
                         
                         PythiaParameters = cms.PSet(
                                                     pythia8CommonSettingsBlock,
                                                     pythia8CP5SettingsBlock,
                                                     processParameters = cms.vstring(
                                                                                     'SUSY:all = off',
                                                                                     'SUSY:qqbar2chi+chi- = on',
                                                                                     'SUSY:qqbar2chi+-chi0 = on',
                                                                                     '1000024:isResonance = false',
                                                                                     '1000024:oneChannel = 1 1.0 100 1000022 211',
                                                                                     '1000024:tau0 = %.1f' % CTAU,
                                                                                     'ParticleDecays:tau0Max = %.1f' % (CTAU * 10),
                                                                                     ),
                                                     parameterSets = cms.vstring(
                                                                                 'pythia8CommonSettings',
                                                                                 'pythia8CP5Settings',
                                                                                 'processParameters')
                                                     ),
                         # The following parameters are required by Exotica_HSCP_SIM_cfi:
                         slhaFile = cms.untracked.string('%s' % SLHA_TABLE),   # value not used
                         processFile = cms.untracked.string('SimG4Core/CustomPhysics/data/RhadronProcessList.txt'),
                         useregge = cms.bool(False),
                         hscpFlavor = cms.untracked.string('stau'),
                         massPoint = cms.untracked.int32(MCHI),  # value not used
                         particleFile = cms.untracked.string('DisappTrks/SignalMC/data/geant4/geant4_AMSB_chargino_%sGeV_ctau%scm.slha' % (MCHI, CTAU/10))
                         )

#     Filter setup
# ------------------------
# https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/HepMCCandAlgos/python/genParticles_cfi.py
tmpGenParticles = cms.EDProducer("GenParticleProducer",
                                 saveBarCodes = cms.untracked.bool(True),
                                 src = cms.InputTag("generator","unsmeared"),
                                 abortOnUnknownPDGCode = cms.untracked.bool(False)
                                 )

# https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/RecoJets/Configuration/python/GenJetParticles_cff.py
# https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/RecoMET/Configuration/python/GenMETParticles_cff.py
tmpGenParticlesForJetsNoNu = cms.EDProducer("InputGenJetsParticleSelector",
                                            src = cms.InputTag("tmpGenParticles"),
                                            ignoreParticleIDs = cms.vuint32(
                                                                            1000022, 1000023, 1000024,
                                                                            1000012, 1000014, 1000016,
                                                                            2000012, 2000014, 2000016,
                                                                            1000039, 5100039,
                                                                            4000012, 4000014, 4000016,
                                                                            9900012, 9900014, 9900016,
                                                                            39,12,14,16),
                                            partonicFinalState = cms.bool(False),
                                            excludeResonances = cms.bool(False),
                                            excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
                                            tausAsJets = cms.bool(False)
                                            )

tmpGenMetTrue = cms.EDProducer("GenMETProducer",
                               src = cms.InputTag("tmpGenParticlesForJetsNoNu"),
                               onlyFiducialParticles = cms.bool(False), ## Use only fiducial GenParticles
                               globalThreshold = cms.double(0.0), ## Global Threshold for input objects
                               usePt   = cms.bool(True), ## using Pt instead Et
                               applyFiducialThresholdForFractions   = cms.bool(False),
                               )

genMETfilter1 = cms.EDFilter("CandViewSelector",
                             src = cms.InputTag("tmpGenMetTrue"),
                             cut = cms.string("pt > 200")
                             )

genMETfilter2 = cms.EDFilter("CandViewCountFilter",
                             src = cms.InputTag("genMETfilter1"),
                             minNumber = cms.uint32(1),
                             )

#ProductionFilterSequence = cms.Sequence(generator)

ProductionFilterSequence = cms.Sequence(generator*
                                        tmpGenParticles * tmpGenParticlesForJetsNoNu *
                                        tmpGenMetTrue * genMETfilter1 * genMETfilter2)
