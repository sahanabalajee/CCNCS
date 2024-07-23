#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"
#include "ns3/flow-monitor-module.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "ns3/energy-module.h"
#include "ns3/basic-energy-source-helper.h"
#include "ns3/wifi-radio-energy-model-helper.h"
#include <sys/stat.h>  
#include <sys/types.h> 
#include <errno.h>     
#include <string.h>
#include "ns3/pcap-file-wrapper.h"
   

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("SimpleSimulation");

std::ofstream throughputFile;
std::ofstream outageFile;
std::ofstream metricsFile;
std::ofstream energyFile;
uint64_t lastTotalRx = 0;
Ptr<PacketSink> sink;
bool isOutage = false;
double totalEnergyConsumed = 0.0;
double totalprevEnergyConsumed = 0.0;

double currEnergyConsumed[3] = {0.0,0.0,0.0};
double prevEnergyConsumed[3] = {0.0,0.0,0.0};

std::string GetTrafficType(Ipv4Address sourceAddress) {
    // Define the IP addresses of DDoS nodes
    static Ipv4Address ddosTcpAddress = Ipv4Address("10.0.0.1");
    static Ipv4Address ddosUdpAddress = Ipv4Address("10.0.0.2");

    if (sourceAddress == ddosTcpAddress) {
        return "DDoS_TCP";
    } else if (sourceAddress == ddosUdpAddress) {
        return "DDoS_UDP";
    } else {
        return "Normal";
    }
}

void CreateOutputFolder(const std::string &folderPath) {
    struct stat info;

    // Check if the directory exists
    if (stat(folderPath.c_str(), &info) != 0) {
        // Directory does not exist
        if (errno == ENOENT) {
            // Try to create the directory
            if (mkdir(folderPath.c_str(), 0777) != 0) {
                std::cerr << "Error creating directory " << folderPath << ": " << strerror(errno) << std::endl;
                exit(EXIT_FAILURE);
            }
        } else {
            // Another error occurred
            std::cerr << "Error checking directory " << folderPath << ": " << strerror(errno) << std::endl;
            exit(EXIT_FAILURE);
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        // Path exists but is not a directory
        std::cerr << "Path " << folderPath << " exists but is not a directory" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void CalculateThroughput (Ipv4Address apAddress) {
    Time now = Simulator::Now();
    double cur = (sink->GetTotalRx () - lastTotalRx) * (double) 8 / 1e6;  // Mbps
    throughputFile << now.GetSeconds() << "\t" << cur << std::endl;
    lastTotalRx = sink->GetTotalRx();

    Simulator::Schedule (Seconds (1.0), &CalculateThroughput, apAddress);
}

void CollectMetrics (Ptr<FlowMonitor> monitor, Ptr<Ipv4FlowClassifier> classifier) {
    std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats ();
    for (auto const &flow : stats) {
        Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (flow.first);
        std::string trafficType = GetTrafficType(t.sourceAddress); // Determine the traffic type
        metricsFile << Simulator::Now().GetSeconds() << ","
                    << flow.first << ","
                    << t.sourceAddress << ","
                    << t.destinationAddress << ","
                    << (t.protocol == 6 ? "TCP" : "UDP") << ","
                    << t.sourcePort << ","
                    << t.destinationPort << ","
                    << flow.second.txBytes << ","
                    << flow.second.txBytes * 8.0 / 1e6 << ","
                    << flow.second.rxBytes * 8.0 / 1e6 << ","
                    << flow.second.lostPackets << ","
                    << (flow.second.delaySum.GetSeconds() / flow.second.rxPackets) * 1e3 << ","
                    << (flow.second.jitterSum.GetSeconds() / flow.second.rxPackets) * 1e3 << ","
                    << stats.size() << ","
                    << "0" << ","  // Error rate placeholder
                    << "0" << ","  // Network load placeholder
                    << "0" << ","  // CPU usage placeholder
                    << "0" << ","  // Memory usage placeholder
                    << trafficType << std::endl;  // Traffic type
    }
    Simulator::Schedule (Seconds (1.0), &CollectMetrics, monitor, classifier);
}


void CalculateEnergyConsumption(NodeContainer sensorNodes, NodeContainer ddosNodes, NodeContainer wifiApNode) {
    totalEnergyConsumed = 0.0;
    double currentEnergy = 0.0;

    // Assuming a map to store energy consumption per node
    std::map<Ptr<Node>, double> nodeEnergyConsumed;

    for (size_t i = 0; i < sensorNodes.GetN(); ++i) {
        Ptr<Node> node = sensorNodes.Get(i);
        Ptr<ns3::energy::BasicEnergySource> energySource = DynamicCast<ns3::energy::BasicEnergySource>(node->GetObject<ns3::energy::EnergySourceContainer>()->Get(0));
        if (energySource) {
            double energyConsumed = energySource->GetInitialEnergy() - energySource->GetRemainingEnergy();
            //std::cout << "Energy consumed: " << energyConsumed << " Joules" << std::endl;
            nodeEnergyConsumed[node] = energyConsumed;
            totalEnergyConsumed += energyConsumed;
            currentEnergy+=energyConsumed;
        } else {
            NS_LOG_WARN("Energy source not found for sensor node " << i);
        }
    }

    for (size_t i = 0; i < ddosNodes.GetN(); ++i) {
        Ptr<Node> node = ddosNodes.Get(i);
        Ptr<ns3::energy::BasicEnergySource> energySource = DynamicCast<ns3::energy::BasicEnergySource>(node->GetObject<ns3::energy::EnergySourceContainer>()->Get(0));
        if (energySource) {
            double energyConsumed = energySource->GetInitialEnergy() - energySource->GetRemainingEnergy();
            nodeEnergyConsumed[node] = energyConsumed;
            currentEnergy+=energyConsumed;
        } else {
            NS_LOG_WARN("Energy source not found for DDoS node " << i);
        }
    }

    for (size_t i = 0; i < wifiApNode.GetN(); ++i) {
        Ptr<Node> node = wifiApNode.Get(i);
        Ptr<ns3::energy::BasicEnergySource> energySource = DynamicCast<ns3::energy::BasicEnergySource>(node->GetObject<ns3::energy::EnergySourceContainer>()->Get(0));
        if (energySource) {
            double energyConsumed = energySource->GetInitialEnergy() - energySource->GetRemainingEnergy();
            nodeEnergyConsumed[node] = energyConsumed;
            currentEnergy+=energyConsumed;
        } else {
            NS_LOG_WARN("Energy source not found for WiFi AP node " << i);
        }
    }
    
    //std::cout << "Energy Difference: " << totalEnergyConsumed - totalprevEnergyConsumed << " Joules" << std::endl;
    if (totalEnergyConsumed - totalprevEnergyConsumed > 15.712 ){
    	outageFile << "Outage started at: " << Now().GetSeconds() << " seconds" << std::endl;}
    

    energyFile << Simulator::Now().GetSeconds() << "\t" <<totalEnergyConsumed << "\t" << std::endl;
    totalprevEnergyConsumed = totalEnergyConsumed;
    totalEnergyConsumed = 0.0;

    Simulator::Schedule(Seconds(1.0), &CalculateEnergyConsumption, sensorNodes, ddosNodes, wifiApNode);
}

int main (int argc, char *argv[]) {
    Time simulationTime = Seconds (60);
    std::string phyMode ("DsssRate11Mbps");
    uint32_t numSensorNodes = 10;  // Increase the number of sensor nodes
    std::string outputFolder = "output/";

    CommandLine cmd;
    cmd.AddValue ("numSensorNodes", "Number of sensor nodes", numSensorNodes);
    cmd.Parse (argc,argv);
    
    CreateOutputFolder(outputFolder);
    CreateOutputFolder(outputFolder+ "pcap");

    NodeContainer wifiApNode;
    wifiApNode.Create (1);
    NodeContainer sensorNodes;
    sensorNodes.Create (numSensorNodes);
    NodeContainer ddosNodes;
    ddosNodes.Create (2);  // Two attacker nodes, one for UDP and one for TCP attacks

    MobilityHelper mobility;
    Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator> ();
    positionAlloc->Add (Vector (0.0, 0.0, 0.0));  // AP position
    for (uint32_t i = 0; i < sensorNodes.GetN (); ++i) {
        positionAlloc->Add (Vector (5.0 * i, 0.0, 0.0));  // Sensor nodes positions
    }
    for (uint32_t i = 0; i < ddosNodes.GetN (); ++i) {
        positionAlloc->Add (Vector (5.0 * (i + sensorNodes.GetN()), 0.0, 0.0));  // DDoS attacker nodes positions
    }

    mobility.SetPositionAllocator (positionAlloc);
    mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
    mobility.Install (wifiApNode);
    mobility.Install (sensorNodes);
    mobility.Install (ddosNodes);

    InternetStackHelper stack;
    stack.Install (wifiApNode);
    stack.Install (sensorNodes);
    stack.Install (ddosNodes);

    YansWifiChannelHelper wifiChannel = YansWifiChannelHelper::Default ();
    YansWifiPhyHelper wifiPhy = YansWifiPhyHelper ();
    wifiPhy.SetChannel (wifiChannel.Create ());

    WifiHelper wifi;
    wifi.SetStandard (WIFI_STANDARD_80211b);
    wifi.SetRemoteStationManager ("ns3::AarfWifiManager");

    WifiMacHelper wifiMac;
    Ssid ssid = Ssid ("network");
    wifiMac.SetType ("ns3::ApWifiMac", "Ssid", SsidValue (ssid));
    NetDeviceContainer apDevices;
    apDevices = wifi.Install (wifiPhy, wifiMac, wifiApNode);

    wifiMac.SetType ("ns3::StaWifiMac", "Ssid", SsidValue (ssid), "ActiveProbing", BooleanValue (false));
    NetDeviceContainer sensorDevices;
    sensorDevices = wifi.Install (wifiPhy, wifiMac, sensorNodes);
    NetDeviceContainer ddosDevices;
    ddosDevices = wifi.Install (wifiPhy, wifiMac, ddosNodes);

    Ipv4AddressHelper ipv4;
    ipv4.SetBase ("10.0.0.0", "255.255.255.0");
    Ipv4InterfaceContainer apInterface = ipv4.Assign (apDevices);
    Ipv4InterfaceContainer sensorInterface = ipv4.Assign (sensorDevices);
    Ipv4InterfaceContainer ddosInterface = ipv4.Assign (ddosDevices);

    uint16_t port = 9;
    Address apLocalAddress (InetSocketAddress (Ipv4Address::GetAny (), port));
    PacketSinkHelper packetSinkHelper ("ns3::TcpSocketFactory", apLocalAddress);
    ApplicationContainer apApp = packetSinkHelper.Install (wifiApNode.Get (0));
    apApp.Start (Seconds (0.0));
    apApp.Stop (simulationTime);

    sink = StaticCast<PacketSink> (apApp.Get (0));

    OnOffHelper onOffHelper ("ns3::TcpSocketFactory", Ipv4Address::GetAny ());
    onOffHelper.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
    onOffHelper.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
    onOffHelper.SetAttribute ("DataRate", DataRateValue (DataRate ("5Mbps")));
    onOffHelper.SetAttribute ("PacketSize", UintegerValue (1024));

    for (uint32_t i = 0; i < sensorNodes.GetN (); ++i) {
        AddressValue remoteAddress (InetSocketAddress (apInterface.GetAddress (0), port));
        onOffHelper.SetAttribute ("Remote", remoteAddress);
        ApplicationContainer tempApp = onOffHelper.Install (sensorNodes.Get (i));
        tempApp.Start (Seconds (1.0));
        tempApp.Stop (simulationTime);
    }

    // Configure DDoS attackers
    OnOffHelper ddosTcpHelper ("ns3::TcpSocketFactory", Address ());
    ddosTcpHelper.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
    ddosTcpHelper.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
    ddosTcpHelper.SetAttribute ("DataRate", DataRateValue (DataRate ("50Mbps")));  // Higher rate for DDoS
    ddosTcpHelper.SetAttribute ("PacketSize", UintegerValue (1024));

    OnOffHelper ddosUdpHelper ("ns3::UdpSocketFactory", Address ("10.0.0.255"));
    ddosUdpHelper.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
    ddosUdpHelper.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
    ddosUdpHelper.SetAttribute ("DataRate", DataRateValue (DataRate ("50Mbps")));  // Higher rate for DDoS
    ddosUdpHelper.SetAttribute ("PacketSize", UintegerValue (1024));

    Ptr<UniformRandomVariable> randomVar = CreateObject<UniformRandomVariable> ();
    randomVar->SetAttribute ("Min", DoubleValue (1.0));
    randomVar->SetAttribute ("Max", DoubleValue (10.0));
   

    for (uint32_t i = 0; i < ddosNodes.GetN (); ++i) {
        AddressValue remoteAddress (InetSocketAddress (apInterface.GetAddress (0), port));
        if (i % 2 == 0) {
            ddosTcpHelper.SetAttribute ("Remote", remoteAddress);
            ApplicationContainer ddosTcpApp = ddosTcpHelper.Install (ddosNodes.Get (i));
            ddosTcpApp.Start (Seconds (randomVar->GetValue ()));
            ddosTcpApp.Stop (simulationTime);
        } else {
            ddosUdpHelper.SetAttribute ("Remote", remoteAddress);
            ApplicationContainer ddosUdpApp = ddosUdpHelper.Install (ddosNodes.Get (i));
            ddosUdpApp.Start (Seconds (randomVar->GetValue ()));
            ddosUdpApp.Stop (simulationTime);
        }
    }
    
    
       
    BasicEnergySourceHelper basicSourceHelper;
    basicSourceHelper.Set("BasicEnergySourceInitialEnergyJ", DoubleValue(1000.0));
    basicSourceHelper.Set("BasicEnergySupplyVoltageV", DoubleValue(12.0));
    
    WifiRadioEnergyModelHelper radioEnergyHelper;
    
    radioEnergyHelper.Set("TxCurrentA",DoubleValue(0.017));
    radioEnergyHelper.Set("RxCurrentA",DoubleValue(0.1));
    radioEnergyHelper.Set("IdleCurrentA",DoubleValue(0.273));
    radioEnergyHelper.Set("SleepCurrentA",DoubleValue(0.033));
    
    ns3::energy::EnergySourceContainer sensorSources = basicSourceHelper.Install(sensorNodes);
    ns3::energy::DeviceEnergyModelContainer sensorModels = radioEnergyHelper.Install(sensorDevices, sensorSources);

    ns3::energy::EnergySourceContainer ddosSources = basicSourceHelper.Install(ddosNodes);
    ns3::energy::DeviceEnergyModelContainer ddosModels = radioEnergyHelper.Install(ddosDevices, ddosSources);

    ns3::energy::EnergySourceContainer apSource = basicSourceHelper.Install(wifiApNode);
    ns3::energy::DeviceEnergyModelContainer apModel = radioEnergyHelper.Install(apDevices, apSource);
    
    Simulator::Schedule(Seconds(1.0), &CalculateEnergyConsumption, sensorNodes,ddosNodes,wifiApNode);
	
    FlowMonitorHelper flowmonHelper;
    Ptr<FlowMonitor> monitor = flowmonHelper.InstallAll ();
    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmonHelper.GetClassifier ());

    throughputFile.open (outputFolder +"throughput.txt");
    outageFile.open (outputFolder +"outage.txt");
    metricsFile.open (outputFolder +"network_metrics.txt");
    energyFile.open (outputFolder +"energy.txt");
    metricsFile << "Timestamp,FlowID,SourceIP,DestinationIP,Protocol,SourcePort,DestinationPort,PacketSize,DataRate,Throughput,PacketLoss,Delay,Jitter,FlowCount,ErrorRate,NetworkLoad,CPUUsage,MemoryUsage,TrafficType" << std::endl;

    Simulator::Schedule (Seconds (1.0), &CalculateThroughput, apInterface.GetAddress (0));
    Simulator::Schedule (Seconds (1.0), &CollectMetrics, monitor, classifier);
    
    //Pcap
	wifiPhy.SetPcapDataLinkType (WifiPhyHelper::DLT_IEEE802_11_RADIO);
	wifiPhy.EnablePcap (outputFolder + "pcap/apDevice", apDevices);
	wifiPhy.EnablePcap (outputFolder + "pcap/sensorDevice", sensorDevices);
	wifiPhy.EnablePcap (outputFolder + "pcap/ddosDevice", ddosDevices);
    

    Simulator::Stop (simulationTime);
    Simulator::Run ();
    Simulator::Destroy ();

    throughputFile.close ();
    outageFile.close ();
    energyFile.close ();
    metricsFile.close ();
    
    monitor->SerializeToXmlFile(outputFolder + "flowmonitor.xml", true, true);
    Simulator::Destroy ();

    return 0;
}
