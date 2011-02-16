#include "Fireworks/Core/interface/CmsShowMainBase.h"

#include "Fireworks/Core/interface/ActionsList.h"
#include "Fireworks/Core/interface/CSGAction.h"
#include "Fireworks/Core/interface/CSGContinuousAction.h"
#include "Fireworks/Core/interface/CmsShowMainFrame.h"
#include "Fireworks/Core/interface/CmsShowSearchFiles.h"
#include "Fireworks/Core/interface/Context.h"
#include "Fireworks/Core/interface/FWColorManager.h"
#include "Fireworks/Core/interface/FWConfigurationManager.h"
#include "Fireworks/Core/interface/FWEveViewManager.h"
#include "Fireworks/Core/interface/FWEventItemsManager.h"
#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/interface/FWJobMetadataManager.h"
#include "Fireworks/Core/interface/FWMagField.h"
#include "Fireworks/Core/interface/FWModelChangeManager.h"
#include "Fireworks/Core/interface/FWNavigatorBase.h"
#include "Fireworks/Core/interface/FWSelectionManager.h"
#include "Fireworks/Core/interface/FWTableViewManager.h"
#include "Fireworks/Core/interface/FWTriggerTableViewManager.h"
#include "Fireworks/Core/interface/FWViewManagerManager.h"
#include "Fireworks/Core/src/CmsShowTaskExecutor.h"
#include "Fireworks/Core/src/FWColorSelect.h"
#include "Fireworks/Core/src/SimpleSAXParser.h"
#include "Fireworks/Core/interface/CmsShowCommon.h"

#include "Fireworks/Core/interface/fwLog.h"

#include "TGLWidget.h"
#include "TGMsgBox.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TTimer.h"
#include "TEveManager.h"

#include <boost/bind.hpp>

CmsShowMainBase::CmsShowMainBase()
   : 
     m_changeManager(new FWModelChangeManager),
     m_colorManager( new FWColorManager(m_changeManager.get())),
     m_configurationManager(new FWConfigurationManager),
     m_eiManager(new FWEventItemsManager(m_changeManager.get())),
     m_guiManager(0),
     m_selectionManager(new FWSelectionManager(m_changeManager.get())),
     m_startupTasks(new CmsShowTaskExecutor),
     m_viewManager(new FWViewManagerManager(m_changeManager.get(), m_colorManager.get())),
     m_autoLoadTimer(new SignalTimer()),
     m_navigatorPtr(0),
     m_metadataManagerPtr(0),
     m_contextPtr(0),
     m_autoLoadTimerRunning(kFALSE),
     m_forward(true),
     m_isPlaying(false),
     m_loop(false),
     m_playDelay(3.f)
{
}

CmsShowMainBase::~CmsShowMainBase()
{
   //avoids a seg fault from eve which happens if eve is terminated after the GUI is gone
   m_selectionManager->clearSelection();
}

void
CmsShowMainBase::setupActions()
{ 
   m_guiManager->writeToPresentConfigurationFile_.connect(sigc::mem_fun(*this, &CmsShowMainBase::writeToCurrentConfigFile));

   // init TGSlider state before signals are connected
   m_guiManager->setDelayBetweenEvents(m_playDelay);

   m_navigatorPtr->newEvent_.connect(boost::bind(&FWGUIManager::loadEvent, guiManager()));
   if (m_guiManager->getAction(cmsshow::sNextEvent) != 0)
      m_guiManager->getAction(cmsshow::sNextEvent)->activated.connect(sigc::mem_fun(*this, &CmsShowMainBase::doNextEvent));
   if (m_guiManager->getAction(cmsshow::sPreviousEvent) != 0)
      m_guiManager->getAction(cmsshow::sPreviousEvent)->activated.connect(sigc::mem_fun(*this, &CmsShowMainBase::doPreviousEvent));
   if (m_guiManager->getAction(cmsshow::sGotoFirstEvent) != 0)
      m_guiManager->getAction(cmsshow::sGotoFirstEvent)->activated.connect(sigc::mem_fun(*this, &CmsShowMainBase::doFirstEvent));
   if (m_guiManager->getAction(cmsshow::sGotoLastEvent) != 0)
      m_guiManager->getAction(cmsshow::sGotoLastEvent)->activated.connect(sigc::mem_fun(*this, &CmsShowMainBase::doLastEvent));
   if (m_guiManager->getAction(cmsshow::sQuit) != 0) 
      m_guiManager->getAction(cmsshow::sQuit)->activated.connect(sigc::mem_fun(*this, &CmsShowMainBase::quit));
 
   m_guiManager->changedEventId_.connect(boost::bind(&CmsShowMainBase::goToRunEvent,this,_1,_2,_3));
   m_guiManager->playEventsAction()->started_.connect(sigc::mem_fun(*this, &CmsShowMainBase::playForward));
   m_guiManager->playEventsBackwardsAction()->started_.connect(sigc::mem_fun(*this,&CmsShowMainBase::playBackward));
   m_guiManager->loopAction()->started_.connect(sigc::mem_fun(*this,&CmsShowMainBase::setPlayLoopImp));
   m_guiManager->loopAction()->stopped_.connect(sigc::mem_fun(*this,&CmsShowMainBase::unsetPlayLoopImp));
   m_guiManager->changedDelayBetweenEvents_.connect(boost::bind(&CmsShowMainBase::setPlayDelay,this,_1));
   m_guiManager->playEventsAction()->stopped_.connect(sigc::mem_fun(*this,&CmsShowMainBase::stopPlaying));
   m_guiManager->playEventsBackwardsAction()->stopped_.connect(sigc::mem_fun(*this,&CmsShowMainBase::stopPlaying));
 
   m_autoLoadTimer->timeout_.connect(boost::bind(&CmsShowMainBase::autoLoadNewEvent, this));
}

void
CmsShowMainBase::setupViewManagers()
{
   guiManager()->updateStatus("Setting up view manager...");

   boost::shared_ptr<FWViewManagerBase> eveViewManager(new FWEveViewManager(guiManager()));
   eveViewManager->setContext(m_contextPtr);
   viewManager()->add(eveViewManager);

   boost::shared_ptr<FWTableViewManager> tableViewManager(new FWTableViewManager(guiManager()));
   configurationManager()->add(std::string("Tables"), tableViewManager.get());
   viewManager()->add(tableViewManager);
   eiManager()->goingToClearItems_.connect(boost::bind(&FWTableViewManager::removeAllItems, tableViewManager.get()));

   boost::shared_ptr<FWTriggerTableViewManager> triggerTableViewManager(new FWTriggerTableViewManager(guiManager()));
   configurationManager()->add(std::string("TriggerTables"), triggerTableViewManager.get());
   configurationManager()->add(std::string("L1TriggerTables"), triggerTableViewManager.get()); // AMT: added for backward compatibilty
   triggerTableViewManager->setContext(m_contextPtr);
   viewManager()->add(triggerTableViewManager);

    
   // Unfortunately, due to the plugin mechanism, we need to delay
   // until here the creation of the FWJobMetadataManager, because
   // otherwise the supportedTypesAndRepresentations map is empty.
   // FIXME: should we have a signal for whenever the above mentioned map
   //        changes? Can that actually happer (maybe if we add support
   //        for loading plugins on the fly??).
   m_metadataManagerPtr->initReps(viewManager()->supportedTypesAndRepresentations());
}

void
CmsShowMainBase::doFirstEvent()
{
   m_navigatorPtr->firstEvent();
   checkPosition();
   draw();
}

void
CmsShowMainBase::doNextEvent()
{
   m_navigatorPtr->nextEvent();
   checkPosition();
   draw();
}

void
CmsShowMainBase::doPreviousEvent()
{
   m_navigatorPtr->previousEvent();
   checkPosition();
   draw();
}
void
CmsShowMainBase::doLastEvent()
{
   m_navigatorPtr->lastEvent();
   checkPosition();
   draw();
}

void
CmsShowMainBase::goToRunEvent(edm::RunNumber_t run, edm::LuminosityBlockNumber_t lumi, edm::EventNumber_t event)
{
   m_navigatorPtr->goToRunEvent(run, lumi, event);
   checkPosition();
   draw();
}


void
CmsShowMainBase::draw()
{
   m_guiManager->updateStatus("loading event ...");

   if (m_contextPtr->getField()->getSource() != FWMagField::kUser)
      m_contextPtr->getField()->checkFiledInfo(m_navigatorPtr->getCurrentEvent());

   TStopwatch sw;
   m_viewManager->eventBegin();
   m_eiManager->newEvent(m_navigatorPtr->getCurrentEvent());
   m_viewManager->eventEnd();
   sw.Stop();
   fwLog(fwlog::kDebug) << "CmsShowMainBase::draw CpuTime " << sw.CpuTime()
                        <<" RealTime " << sw.RealTime() << std::endl;

   if (!m_autoSaveAllViewsFormat.empty())
   {
      m_guiManager->updateStatus("auto saving images ...");
      m_guiManager->exportAllViews(m_autoSaveAllViewsFormat);
   }

   m_guiManager->clearStatus();
}

void
CmsShowMainBase::setup(FWNavigatorBase *navigator,
                       fireworks::Context *context,
                       FWJobMetadataManager *metadataManager)
{
   m_navigatorPtr = navigator;
   m_contextPtr = context;
   m_metadataManagerPtr = metadataManager;
   
   m_colorManager->initialize();
   m_contextPtr->initEveElements();
   m_guiManager.reset(new FWGUIManager(m_contextPtr, m_viewManager.get(), m_navigatorPtr));

   m_eiManager->newItem_.connect(boost::bind(&FWModelChangeManager::newItemSlot,
                                             m_changeManager.get(), _1) );
   
   m_eiManager->newItem_.connect(boost::bind(&FWViewManagerManager::registerEventItem,
                                             m_viewManager.get(), _1));
   m_configurationManager->add("EventItems",m_eiManager.get());
   m_configurationManager->add("GUI",m_guiManager.get());
   m_configurationManager->add("EventNavigator", m_navigatorPtr);
   m_configurationManager->add("CommonPreferences", m_contextPtr->commonPrefs()); // must be after GUIManager in alphabetical order

   m_guiManager->writeToConfigurationFile_.connect(boost::bind(&CmsShowMainBase::writeToConfigFile,
                                                                this,_1));

   m_guiManager->loadFromConfigurationFile_.connect(boost::bind(&CmsShowMainBase::reloadConfiguration,
                                                                this, _1));
   std::string macPath(gSystem->Getenv("CMSSW_BASE"));
   macPath += "/src/Fireworks/Core/macros";
   const char* base = gSystem->Getenv("CMSSW_RELEASE_BASE");
   if(0!=base) {
      macPath+=":";
      macPath +=base;
      macPath +="/src/Fireworks/Core/macros";
   }
   gROOT->SetMacroPath((std::string("./:")+macPath).c_str());
   
   m_startupTasks->tasksCompleted_.connect(boost::bind(&FWGUIManager::clearStatus,
                                                       m_guiManager.get()) );
}

void
CmsShowMainBase::writeToConfigFile(const std::string &name)
{
   m_configFileName = name;
   m_configurationManager->writeToFile(m_configFileName);
}

void
CmsShowMainBase::writeToCurrentConfigFile()
{
   m_configurationManager->writeToFile(m_configFileName);
}

void
CmsShowMainBase::reloadConfiguration(const std::string &config)
{
   if (config.empty())
      return;

   m_configFileName = config;

   std::string msg = "Reloading configuration "
      + config + "...";
   fwLog(fwlog::kDebug) << msg << std::endl;
   m_guiManager->updateStatus(msg.c_str());
   m_guiManager->subviewDestroyAll();
   m_eiManager->clearItems();
   m_configFileName = config;
   try
   {
      gEve->DisableRedraw();
      m_configurationManager->readFromFile(config);
      gEve->EnableRedraw();
   }
   catch (std::runtime_error &e)
   {
      Int_t chosen;
      new TGMsgBox(gClient->GetDefaultRoot(),
                   gClient->GetDefaultRoot(),
                   "Bad configuration",
                   ("Configuration " + config + " cannot be parsed.").c_str(),
                   kMBIconExclamation,
                   kMBCancel,
                   &chosen);
   }
   catch (SimpleSAXParser::ParserError &e)
   {
      Int_t chosen;
      new TGMsgBox(gClient->GetDefaultRoot(),
                   gClient->GetDefaultRoot(),
                   "Bad configuration",
                   ("Configuration " + config + " cannot be parsed.").c_str(),
                   kMBIconExclamation,
                   kMBCancel,
                   &chosen);
   }
   m_guiManager->updateStatus("");
}

void
CmsShowMainBase::setupAutoLoad(float x)
{
   m_playDelay = x;
   m_guiManager->setDelayBetweenEvents(m_playDelay);
   if (!m_guiManager->playEventsAction()->isEnabled())
      m_guiManager->playEventsAction()->enable();

   m_guiManager->playEventsAction()->switchMode();
}

void
CmsShowMainBase::startAutoLoadTimer()
{
   m_autoLoadTimer->SetTime((Long_t)(m_playDelay*1000));
   m_autoLoadTimer->Reset();
   m_autoLoadTimer->TurnOn();
   m_autoLoadTimerRunning = kTRUE;
}

void 
CmsShowMainBase::stopAutoLoadTimer()
{
   m_autoLoadTimer->TurnOff();
   m_autoLoadTimerRunning = kFALSE;
}

void
CmsShowMainBase::setupConfiguration()
{
   m_guiManager->updateStatus("Setting up configuration...");
   if(m_configFileName.empty() ) {
      fwLog(fwlog::kInfo) << "no configuration is loaded." << std::endl;
      m_guiManager->getMainFrame()->MapSubwindows();
      m_guiManager->getMainFrame()->Layout();
      m_guiManager->getMainFrame()->MapRaised();
      m_configFileName = "newconfig.fwc";
      m_guiManager->createView("Rho Phi"); 
      m_guiManager->createView("Rho Z"); 
   }
   else {
      char* whereConfig = gSystem->Which(TROOT::GetMacroPath(), m_configFileName.c_str(), kReadPermission);
      if(0==whereConfig) {
         fwLog(fwlog::kInfo) <<"unable to load configuration file '"<<m_configFileName<<"' will load default instead."<<std::endl;
         whereConfig = gSystem->Which(TROOT::GetMacroPath(), "default.fwc", kReadPermission);
         assert(whereConfig && "Default configuration cannot be found. Malformed Fireworks installation?");
      }
      m_configFileName = whereConfig;

      delete [] whereConfig;
      try
      {
         gEve->DisableRedraw();
         m_configurationManager->readFromFile(m_configFileName);
         gEve->EnableRedraw();
      }
      catch (std::runtime_error &e)
      {
         fwLog(fwlog::kError) <<"Unable to load configuration file '" 
                              << m_configFileName 
                              << "' which was specified on command line. Quitting." 
                              << std::endl;
         exit(1);
      }
      catch (SimpleSAXParser::ParserError &e)
      {
         fwLog(fwlog::kError) <<"Unable to load configuration file '" 
                              << m_configFileName 
                              << "' which was specified on command line. Quitting." 
                              << std::endl;
         exit(1);
      }
   }
   /* 
      if(not m_configFileName.empty() ) {
      //when the program quits we will want to save the configuration automatically
      m_guiManager->goingToQuit_.connect(
      boost::bind(&FWConfigurationManager::writeToFile,
      m_configurationManager.get(),
      m_configFileName));
      }
   */
}


void
CmsShowMainBase::setPlayDelay(Float_t val)
{
   m_playDelay = val;
}

void
CmsShowMainBase::setupDebugSupport()
{
   m_guiManager->updateStatus("Setting up Eve debug window...");
   m_guiManager->openEveBrowserForDebugging();
}

void
CmsShowMainBase::setPlayLoop()
{
   if(!m_loop) {
      setPlayLoopImp();
      m_guiManager->loopAction()->activated();
   }
}

void
CmsShowMainBase::unsetPlayLoop()
{
   if(m_loop) {
      unsetPlayLoopImp();
      m_guiManager->loopAction()->stop();
   }
}

void
CmsShowMainBase::setPlayLoopImp()
{
   m_loop = true;
}

void
CmsShowMainBase::unsetPlayLoopImp()
{
   m_loop = false;
}

void 
CmsShowMainBase::registerPhysicsObject(const FWPhysicsObjectDesc&iItem)
{
   m_eiManager->add(iItem);
}

void
CmsShowMainBase::playForward()
{
   m_forward = true;
   m_isPlaying = true;
   guiManager()->enableActions(kFALSE);
   startAutoLoadTimer();
}

void
CmsShowMainBase::playBackward()
{
   m_forward = false;
   m_isPlaying = true;
   guiManager()->enableActions(kFALSE);
   startAutoLoadTimer();
}

void
CmsShowMainBase::loadGeometry()
{   // prepare geometry service
   // ATTN: this should be made configurable
   try 
   {
      guiManager()->updateStatus("Loading geometry...");
      m_geom.loadMap(m_geometryFilename.c_str());
      m_contextPtr->setGeom(&m_geom);
   }
   catch (const std::runtime_error& iException)
   {
      fwLog(fwlog::kError) << "CmsShowMain::loadGeometry() caught exception: \n"
                           << m_geometryFilename << " "
                           << iException.what() << std::endl;
      exit(0);
   }
}
