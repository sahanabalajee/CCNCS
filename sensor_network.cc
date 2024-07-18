#include "ns3/command-line.h"
#include "ns3/config.h"
#include "ns3/internet-stack-helper.h"
#include "ns3/ipv4-address-helper.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/log.h"
#include "ns3/mobility-helper.h"
#include "ns3/mobility-model.h"
#include "ns3/on-off-helper.h"
#include "ns3/packet-sink-helper.h"
#include "ns3/packet-sink.h"
#include "ns3/ssid.h"
#include "ns3/string.h"
#include "ns3/tcp-westwood-plus.h"
#include "ns3/yans-wifi-channel.h"
#include "ns3/yans-wifi-helper.h"
#include "ns3/flow-monitor.h"
#include "ns3/flow-monitor-helper.h"
#include "ns3/internet-module.h"
#include "ns3/ipv4-flow-classifier.h"
#include "ns3/netanim-module.h"
#include "ns3/constant-velocity-mobility-model.h"
#include "ns3/energy-module.h"
#include "ns3/wifi-radio-energy-model-helper.h"



#include <fstream>

NS_LOG_COMPONENT_DEFINE("proj");

using namespace ns3;

Ptr<PacketSink> sink;
Ptr<PacketSink> sink2;
uint64_t lastTotalRx = 0;
std::ofstream throughputFile;
double totalEnergyConsumed = 0.0;
std::vector<double> nodeEnergyConsumed;


void
CalculateThroughput()
{
    Time now = Simulator::Now(); /* Return the simulator's virtual time. */
    double cur = (sink->GetTotalRx() - lastTotalRx) * 8.0 /1e6; /* Convert Application RX Packets to MBits. */
                
    throughputFile << now.GetSeconds() << "\t" << cur << std::endl;
    

    lastTotalRx = sink->GetTotalRx();
    Simulator::Schedule(MilliSeconds(1000), &CalculateThroughput);
}
void CalculateEnergyConsumption(NodeContainer temperatureSensorNodes, NodeContainer humiditySensorNodes,NodeContainer pressureSensorNodes,NodeContainer soundSensorNodes) {
    totalEnergyConsumed = 0.0;

    for (size_t i = 0; i < nodeEnergyConsumed.size(); ++i) {
        Ptr<ns3::energy::BasicEnergySource> energySource = DynamicCast<ns3::energy::BasicEnergySource>(temperatureSensorNodes.Get(i)->GetObject<ns3::energy::EnergySourceContainer>()->Get(0));
        if (energySource) {
            nodeEnergyConsumed[i] = energySource->GetInitialEnergy() - energySource->GetRemainingEnergy();
            totalEnergyConsumed += nodeEnergyConsumed[i];
        } else {
            NS_LOG_WARN("Energy source not found for temperature sensor node " << i);
        }
    }

    for (size_t i = 0; i < humiditySensorNodes.GetN(); ++i) {
        Ptr<ns3::energy::BasicEnergySource> energySource = DynamicCast<ns3::energy::BasicEnergySource>(humiditySensorNodes.Get(i)->GetObject<ns3::energy::EnergySourceContainer>()->Get(0));
        if (energySource) {
            double energyConsumed = energySource->GetInitialEnergy() - energySource->GetRemainingEnergy();
            totalEnergyConsumed += energyConsumed;
        } else {
            NS_LOG_WARN("Energy source not found for humidity sensor node " << i);
        }
    }
    
        for (size_t i = 0; i < humiditySensorNodes.GetN(); ++i) {
        Ptr<ns3::energy::BasicEnergySource> energySource = DynamicCast<ns3::energy::BasicEnergySource>(pressureSensorNodes.Get(i)->GetObject<ns3::energy::EnergySourceContainer>()->Get(0));
        if (energySource) {
            double energyConsumed = energySource->GetInitialEnergy() - energySource->GetRemainingEnergy();
            totalEnergyConsumed += energyConsumed;
        } else {
            NS_LOG_WARN("Energy source not found for humidity sensor node " << i);
        }
    }
    
        for (size_t i = 0; i < humiditySensorNodes.GetN(); ++i) {
        Ptr<ns3::energy::BasicEnergySource> energySource = DynamicCast<ns3::energy::BasicEnergySource>(soundSensorNodes.Get(i)->GetObject<ns3::energy::EnergySourceContainer>()->Get(0));
        if (energySource) {
            double energyConsumed = energySource->GetInitialEnergy() - energySource->GetRemainingEnergy();
            totalEnergyConsumed += energyConsumed;
        } else {
            NS_LOG_WARN("Energy source not found for humidity sensor node " << i);
        }
    }
    
    

    Simulator::Schedule(Seconds(1.0), &CalculateEnergyConsumption, temperatureSensorNodes, humiditySensorNodes,pressureSensorNodes,soundSensorNodes);
}



int
main(int argc, char *argv[]){
    std::string tcpVariant{"TcpCubic"}; /* TCP variant type. */
    std::string phyRate{"HtMcs7"};        /* Physical layer bitrate. */
    Time simulationTime{"10s"};           /* Simulation time. */

    /* Command line argument parser setup. */
    CommandLine cmd(__FILE__);
    cmd.AddValue("tcpVariant",
                 "Transport protocol to use: TcpNewReno, "
                 "TcpHybla, TcpHighSpeed, TcpHtcp, TcpVegas, TcpScalable, TcpVeno, "
                 "TcpBic, TcpYeah, TcpIllinois, TcpWestwood, TcpWestwoodPlus, TcpLedbat ",
                 tcpVariant);
    cmd.AddValue("phyRate", "Physical layer bitrate", phyRate);
    cmd.AddValue("simulationTime", "Simulation time in seconds", simulationTime);
    cmd.Parse(argc, argv);
    std::string tcpName = tcpVariant;

    tcpVariant = std::string("ns3::") + tcpVariant;
    // Select TCP variant
    TypeId tcpTid;
    NS_ABORT_MSG_UNLESS(TypeId::LookupByNameFailSafe(tcpVariant, &tcpTid),
                        "TypeId " << tcpVariant << " not found");
    Config::SetDefault("ns3::TcpL4Protocol::SocketType",
                       TypeIdValue(TypeId::LookupByName(tcpVariant)));
                       

    /* Configure TCP Options */
   // Config::SetDefault("ns3::TcpSocket::SegmentSize", UintegerValue(payloadSize));
    
    WifiMacHelper wifiMac;
    WifiHelper wifiHelper;
    wifiHelper.SetStandard(WIFI_STANDARD_80211n);

    /* Set up Legacy Channel */
    YansWifiChannelHelper wifiChannel;
    wifiChannel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
    wifiChannel.AddPropagationLoss("ns3::FriisPropagationLossModel", "Frequency", DoubleValue(5e9));

    /* Setup Physical Layer */
    YansWifiPhyHelper wifiPhy;
    wifiPhy.SetChannel(wifiChannel.Create());
    wifiPhy.SetErrorRateModel("ns3::YansErrorRateModel");
    wifiHelper.SetRemoteStationManager("ns3::ConstantRateWifiManager",
                                       "DataMode",
                                       StringValue(phyRate),
                                       "ControlMode",
                                       StringValue("HtMcs0"));
                                       

    NodeContainer apWifiNode;
    int number_of_ap = 4;
    apWifiNode.Create(number_of_ap);
    int number_of_sensors = 21;
    
    NodeContainer temperatureSensorNodes;
    temperatureSensorNodes.Create(number_of_sensors);
    
    NodeContainer humiditySensorNodes;
    humiditySensorNodes.Create(10);
    
    NodeContainer soundSensorNodes;
    soundSensorNodes.Create(10);
    
    NodeContainer pressureSensorNodes;
    pressureSensorNodes.Create(10);
    
    NodeContainer attackerNode;
    attackerNode.Create(1);
    
    Ssid ssid = Ssid("sensor_network");
    wifiMac.SetType("ns3::ApWifiMac", "Ssid", SsidValue(ssid));
    NetDeviceContainer apDevice;
    apDevice = wifiHelper.Install(wifiPhy, wifiMac, apWifiNode);
    
    wifiMac.SetType("ns3::StaWifiMac","Ssid",SsidValue(ssid));
    
    NetDeviceContainer temperatureSensorDevices;
    temperatureSensorDevices = wifiHelper.Install(wifiPhy,wifiMac,temperatureSensorNodes);
    
    
    NetDeviceContainer attackerDevice;
    attackerDevice = wifiHelper.Install(wifiPhy,wifiMac,attackerNode);
    
    NetDeviceContainer humiditySensorDevices;
    humiditySensorDevices = wifiHelper.Install(wifiPhy,wifiMac,humiditySensorNodes);
    
    NetDeviceContainer soundSensorDevices;
    soundSensorDevices = wifiHelper.Install(wifiPhy,wifiMac,soundSensorNodes);
    
    NetDeviceContainer pressureSensorDevices;
    pressureSensorDevices = wifiHelper.Install(wifiPhy,wifiMac,pressureSensorNodes);
    
    
    
    
    InternetStackHelper stack;
    stack.Install(apWifiNode);
    stack.Install(temperatureSensorNodes);
    stack.Install(attackerNode);
    stack.Install(humiditySensorNodes);
    stack.Install(soundSensorNodes);
    stack.Install(pressureSensorNodes);
    
    Ipv4AddressHelper address;
    address.SetBase("192.168.0.0", "255.255.0.0");
    Ipv4InterfaceContainer apInterface;
    apInterface = address.Assign(apDevice);
    Ipv4InterfaceContainer sensorInterface;
    sensorInterface = address.Assign(temperatureSensorDevices);
    Ipv4InterfaceContainer attackerInterface;
    attackerInterface = address.Assign(attackerDevice);
    Ipv4InterfaceContainer humditySensorInterface;
    humditySensorInterface = address.Assign(humiditySensorDevices);
    Ipv4InterfaceContainer soundSensorInterface;
    soundSensorInterface = address.Assign(soundSensorDevices);
    Ipv4InterfaceContainer pressureSensorInterface;
    pressureSensorInterface = address.Assign(pressureSensorDevices);
    
    
    Ptr<UniformRandomVariable> x = CreateObject<UniformRandomVariable>();
    x->SetAttribute ("Min", DoubleValue (112.5));
	x->SetAttribute ("Max", DoubleValue (250.0));
	
	Ptr<UniformRandomVariable> y = CreateObject<UniformRandomVariable>();
    y->SetAttribute ("Min", DoubleValue (0.0));
	y->SetAttribute ("Max", DoubleValue (117.5));
	
	Ptr<UniformRandomVariable> z = CreateObject<UniformRandomVariable>();
    z->SetAttribute ("Min", DoubleValue (0.0));
	z->SetAttribute ("Max", DoubleValue (0.0));
    
    MobilityHelper soundMobility;
    Ptr<ListPositionAllocator> soundAlloc = CreateObject <ListPositionAllocator>();
    
    for(int i = 0; i < 10; i++){
    	double pos_x = x->GetValue();
    	double pos_y = y->GetValue();
    	double pos_z = z->GetValue();
    	
    	soundAlloc->Add(Vector(pos_x,pos_y,pos_z));
    }
    
    
    soundMobility.SetPositionAllocator(soundAlloc);
    soundMobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
	soundMobility.Install(soundSensorNodes);    
    
    MobilityHelper sensorMobility;
    sensorMobility.SetPositionAllocator("ns3::GridPositionAllocator","MinX", DoubleValue(0.0),"MinY",DoubleValue(0.0),"DeltaX",DoubleValue(20),"DeltaY",DoubleValue(20),"LayoutType",StringValue("RowFirst"),"GridWidth",UintegerValue(5));
   
   
   MobilityHelper humidityMobility; 
     Ptr<ListPositionAllocator> humidityPositionAlloc = CreateObject<ListPositionAllocator>();
        
        double radius = 187.5;
        
   for(int i = 0; i < 10; i++){
        	double angle = 2 * M_PI * i/ 10;
        	double x = 62.5 + 50 * cos(angle);
        	double y = radius + 50 * sin(angle);
        	humidityPositionAlloc->Add(Vector(x,y,0.0));
    }
    humidityMobility.SetPositionAllocator(humidityPositionAlloc);
    humidityMobility.Install(humiditySensorNodes);
    
    sensorMobility.Install(temperatureSensorNodes);
    
    
    MobilityHelper apMobility;
    
    Ptr<ListPositionAllocator> apPositionAlloc = CreateObject <ListPositionAllocator>();
    
    apPositionAlloc->Add(Vector(70.0,50.0,0.0));
    apPositionAlloc->Add(Vector(62.5,187.5,0.0));
    apPositionAlloc->Add(Vector(161.3,55.5,0.0));
    apPositionAlloc->Add(Vector(150.0,176.25,0.0));
    
    apMobility.SetPositionAllocator(apPositionAlloc);
    apMobility.Install(apWifiNode); 
    
    
	MobilityHelper attackerMobility;
	Ptr<ListPositionAllocator> attackerPosAlloc = CreateObject <ListPositionAllocator>();
	
	attackerPosAlloc->Add(Vector(30.0,70.0,0.0));
	attackerMobility.SetPositionAllocator(attackerPosAlloc);
	attackerMobility.Install(attackerNode);
    
    
    MobilityHelper pressureMobility;
    Ptr<ListPositionAllocator> pressurePosAlloc = CreateObject <ListPositionAllocator>();
    
    int m = 1;
    
    for(double i = 125; i < 250; i+=25){
    	double y = m*i;
    	pressurePosAlloc->Add(Vector(i,y, 0.0));
    }
    
    for(double i = 125; i <250 ; i+=25){
    	double y = (-m)*i + 350;
    	pressurePosAlloc->Add(Vector(i,y,0.0));
    }
    
    pressureMobility.SetPositionAllocator(pressurePosAlloc);
    pressureMobility.Install(pressureSensorNodes);
    
    Ipv4GlobalRoutingHelper::PopulateRoutingTables();
    
    PacketSinkHelper sinkHelper("ns3::TcpSocketFactory",InetSocketAddress(InetSocketAddress(Ipv4Address::GetAny(), 9)));
    ApplicationContainer sinkApp = sinkHelper.Install(apWifiNode.Get(0));
    ApplicationContainer secondSinkApp = sinkHelper.Install(apWifiNode.Get(1));
    
    sink2 = StaticCast<PacketSink>(secondSinkApp.Get(0));
    sink = StaticCast<PacketSink>(sinkApp.Get(0));
    
    
    OnOffHelper temperaturePktServer("ns3::TcpSocketFactory", (InetSocketAddress(apInterface.GetAddress(0), 9)));
    temperaturePktServer.SetAttribute("PacketSize", UintegerValue(3));
    temperaturePktServer.SetAttribute("OnTime", StringValue("ns3::ConstantRandomVariable[Constant=1]"));
    temperaturePktServer.SetAttribute("OffTime", StringValue("ns3::ConstantRandomVariable[Constant=1]"));
    temperaturePktServer.SetAttribute("DataRate", DataRateValue(DataRate("10Kb/s")));
    ApplicationContainer temperaturePktServerApp = temperaturePktServer.Install(temperatureSensorNodes);
    
    OnOffHelper humidityPktServer("ns3::TcpSocketFactory", (InetSocketAddress(apInterface.GetAddress(1), 9)));
    humidityPktServer.SetAttribute("PacketSize", UintegerValue(2));
    humidityPktServer.SetAttribute("OnTime", StringValue("ns3::ConstantRandomVariable[Constant=1]"));
    humidityPktServer.SetAttribute("OffTime", StringValue("ns3::ConstantRandomVariable[Constant=1]"));
    humidityPktServer.SetAttribute("DataRate", DataRateValue(DataRate("10Kb/s")));
    ApplicationContainer humidityPktServerApp = humidityPktServer.Install(humiditySensorNodes);
    
    OnOffHelper pressurePktServer("ns3::TcpSocketFactory", (InetSocketAddress(apInterface.GetAddress(2), 9)));
    pressurePktServer.SetAttribute("PacketSize", UintegerValue(5));
    pressurePktServer.SetAttribute("OnTime", StringValue("ns3::ConstantRandomVariable[Constant=1]"));
    pressurePktServer.SetAttribute("OffTime", StringValue("ns3::ConstantRandomVariable[Constant=2]"));
    pressurePktServer.SetAttribute("DataRate", DataRateValue(DataRate("20Kb/s")));
    ApplicationContainer pressurePktServerApp = pressurePktServer.Install(pressureSensorNodes);
    
    OnOffHelper soundPktServer("ns3::TcpSocketFactory", (InetSocketAddress(apInterface.GetAddress(2), 9)));
    soundPktServer.SetAttribute("PacketSize", UintegerValue(5));
    soundPktServer.SetAttribute("OnTime", StringValue("ns3::ConstantRandomVariable[Constant=1]"));
    soundPktServer.SetAttribute("OffTime", StringValue("ns3::ConstantRandomVariable[Constant=2]"));
    soundPktServer.SetAttribute("DataRate", DataRateValue(DataRate("20Kb/s")));
    ApplicationContainer soundPktServerApp = soundPktServer.Install(soundSensorNodes);
    
    OnOffHelper tcpDOS("ns3::TcpSocketFactory",(InetSocketAddress(apInterface.GetAddress(0),9)));
    tcpDOS.SetAttribute("PacketSize",UintegerValue(43)); /*Replace ENTER_PAYLOAD_SIZE_IN_BYTES_HERE with the avg payload size for DOS attacks*/
    tcpDOS.SetAttribute("OnTime",StringValue("ns3::ConstantRandomVariable[Constant=100]"));
    tcpDOS.SetAttribute("OffTime",StringValue("ns3::ConstantRandomVariable[Constant=0]"));
    tcpDOS.SetAttribute("DataRate",DataRateValue(DataRate("5Mbps")));
    ApplicationContainer tcpDOSApp = tcpDOS.Install(attackerNode);
    
    
    

    FlowMonitorHelper flowmon;
    Ptr<FlowMonitor> monitor = flowmon.InstallAll();
        
        
    //std :: string throughputFileName = "throughput/throughput_" + tcpName + std:: string {".txt"};
    
    std::string throughputFileName = "throughput.txt";

    
    throughputFile.open(throughputFileName);
    AnimationInterface anim("visual.xml");
   
    for(int i = 0; i < number_of_ap; i++){
    Ptr<MobilityModel> apmob = apWifiNode.Get(i)->GetObject<MobilityModel>();
    double apx = apmob->GetPosition().x;
    double apy = apmob->GetPosition().y;
    anim.SetConstantPosition(apWifiNode.Get(i),apx,apy);
    anim.UpdateNodeColor(apWifiNode.Get(i),0,0,255);
    
    }
    

    
    Ptr<MobilityModel> attackerMob = attackerNode.Get(0)->GetObject<MobilityModel>();
    double ax = attackerMob->GetPosition().x;
    double ay = attackerMob->GetPosition().y;
    
    anim.SetConstantPosition(attackerNode.Get(0),ax,ay);
    anim.UpdateNodeColor(attackerNode.Get(0),0,0,0);
    
    
    sinkApp.Start(Seconds(0.0));
    temperaturePktServerApp.Start(Seconds(1.1));
    humidityPktServerApp.Start(Seconds(1.2));
    soundPktServerApp.Start(Seconds(1.3));
    pressurePktServerApp.Start(Seconds(1.4));
   // tcpDOSApp.Start(Seconds(1.0));
    
    BasicEnergySourceHelper basicSourceHelper;
    basicSourceHelper.Set("BasicEnergySourceInitialEnergyJ", DoubleValue(1000.0));
    basicSourceHelper.Set("BasicEnergySupplyVoltageV", DoubleValue(12.0));
    
    WifiRadioEnergyModelHelper radioEnergyHelper;
    
    radioEnergyHelper.Set("TxCurrentA",DoubleValue(0.017));
    radioEnergyHelper.Set("RxCurrentA",DoubleValue(0.0197));
    radioEnergyHelper.Set("IdleCurrentA",DoubleValue(0.273));
    radioEnergyHelper.Set("SleepCurrentA",DoubleValue(0.033));
    
    ns3::energy::EnergySourceContainer sources = basicSourceHelper.Install(temperatureSensorNodes);
    ns3::energy::DeviceEnergyModelContainer deviceModels = radioEnergyHelper.Install(temperatureSensorDevices, sources);
    
    ns3::energy::EnergySourceContainer srcs_h = basicSourceHelper.Install(humiditySensorNodes);
    ns3::energy::DeviceEnergyModelContainer devModels_h = radioEnergyHelper.Install(humiditySensorDevices,srcs_h);
    
    ns3::energy::EnergySourceContainer srcs_s = basicSourceHelper.Install(soundSensorNodes);
    ns3::energy::DeviceEnergyModelContainer devModels_s = radioEnergyHelper.Install(soundSensorDevices,srcs_s);
    
    ns3::energy::EnergySourceContainer srcs_p = basicSourceHelper.Install(pressureSensorNodes);
    ns3::energy::DeviceEnergyModelContainer devModels_p = radioEnergyHelper.Install(pressureSensorDevices,srcs_p);

    nodeEnergyConsumed.resize(temperatureSensorNodes.GetN(), 0.0);

    Simulator::Schedule(Seconds(1.0), &CalculateEnergyConsumption,temperatureSensorNodes,humiditySensorNodes,pressureSensorNodes,soundSensorNodes);
   
    
    Simulator::Schedule(Seconds(1.1), &CalculateThroughput);
    
    
    
    for(int i = 0; i < number_of_sensors; i++){
    	Ptr<MobilityModel> mob = temperatureSensorNodes.Get(i)->GetObject<MobilityModel>();
    	double x = mob->GetPosition().x;
    	double y = mob->GetPosition().y;
    	anim.SetConstantPosition(temperatureSensorNodes.Get(i),x,y);
        anim.UpdateNodeColor(temperatureSensorNodes.Get(i),0,255,0);
    }
    
    for(int i = 0; i < 10; i++){
    	Ptr<MobilityModel> mob = humiditySensorNodes.Get(i)->GetObject<MobilityModel>();
    	double x = mob->GetPosition().x;
    	double y = mob->GetPosition().y;
    	anim.SetConstantPosition(temperatureSensorNodes.Get(i),x,y);
        anim.UpdateNodeColor(temperatureSensorNodes.Get(i),100,100,10);
    }
    for(int i = 0; i < 10; i++){
    	Ptr<MobilityModel> mob = pressureSensorNodes.Get(i)->GetObject<MobilityModel>();
    	double x = mob->GetPosition().x;
    	double y = mob->GetPosition().y;
    	anim.SetConstantPosition(temperatureSensorNodes.Get(i),x,y);
        anim.UpdateNodeColor(pressureSensorNodes.Get(i),25,50,75);
    }
    
   
    
    
    Simulator::Stop(simulationTime);
    Simulator :: Run();
    
	double averageEnergyConsumption = totalEnergyConsumed / temperatureSensorNodes.GetN();
    std::cout << "Average energy consumption: " << averageEnergyConsumption << " J" << std::endl;
    
    throughputFile.close();
    
    auto averageThroughput =
        (static_cast<double>(sink->GetTotalRx() * 8  ) / simulationTime.GetMicroSeconds());
    
    std::cout << "\nAverage throughput: " << averageThroughput << " Mbit/s" << std::endl;
    
    //Flow monitor code
    monitor->CheckForLostPackets();
    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmon.GetClassifier());
    std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats();

    std::ofstream flowStatsFile;
    
    std::string flowFileName = "flowstats.txt";
     
    flowStatsFile.open(flowFileName);
    
    int total_tx = 0;
    int total_rx = 0;
    
    // Print Flow Monitor statistics for flows terminating at the sink
    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin (); i != stats.end (); ++i)
    {
        Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (i->first);
        if (t.destinationAddress == apInterface.GetAddress (0))
        {
            total_tx+= i->second.txPackets;
            total_rx+= i->second.rxPackets;
            flowStatsFile << i->first << "\t"
                          << t.sourceAddress << "\t"
                          << t.destinationAddress << "\t"
                          << i->second.txBytes << "\t"
                          << i->second.rxBytes << "\t"
                          << i->second.txPackets << "\t"
                          << i->second.rxPackets << "\t"
                          << i->second.lostPackets << "\t"
                          << i->second.rxBytes * 8.0 / (i->second.timeLastRxPacket.GetSeconds() - i->second.timeFirstTxPacket.GetSeconds()) / 1024 / 1024 << std::endl;
       
        }
    }
    

    flowStatsFile.close();
    
    std ::ofstream pktStatsFile;
    
    pktStatsFile.open("packet_stats.txt",std::ios::app);
    
    pktStatsFile << tcpName << "\t" << total_tx << "\t" << total_rx << std::endl;
    
    pktStatsFile.close();
   
    Simulator :: Destroy();
    
    std::string energyStats = "energystats/energy_"+tcpName+".txt";
    
    std::ofstream energyFile;
    energyFile.open(energyStats);
    
    for (size_t i = 0; i < nodeEnergyConsumed.size(); ++i) {
    	double power_consumed = nodeEnergyConsumed[i] / 10;
        energyFile << "Node " << i << ": " << power_consumed << " W" << std::endl;
    }
    energyFile.close();

  
    
    return 0;
}
