# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import geopandas as gpd
from shapely import wkt, geometry
import numpy as np
import time
import csv
import sys
import os



class Node:
    def __init__(self, node_id):
        self.node_id = node_id
        self.x_coord = None
        self.y_coord = None
        self.centroid = None
        self.zone_id = None
        self.geometry = None
        
        self.other_attrs = {}
   

   
class Link:
    def __init__(self, link_id):
        self.org_link_id = None
        self.link_id = link_id
        self.from_node = None
        self.to_node = None
        self.lanes = None
        self.length = 0
        self.dir_flag = 1
        self.geometry = None
        
        self.free_speed = None
        self.capacity = None
        self.link_type = 0
        self.vdf_code = 0
        
        
        self.other_attrs = {}

class Network:
    def __init__(self):
        self.node_dict = {}
        self.link_dict = {}

        self.max_node_id = 0
        self.max_link_id = 0

        self.original_coordinate_type = 'lonlat'

        self.node_other_attrs = []
        self.link_other_attrs = []

    @property
    def number_of_nodes(self):
        return len(self.node_dict)

    @property
    def number_of_links(self):
        return len(self.link_dict)

   

def loadCSVFromSHP(shapfile_path):
    start_time = time.process_time()

    print('Loading shapefile with geometry ...')
    raw_network_csv = gpd.read_file(shapfile_path)
    raw_network_csv.plot()
    
    end_time = time.process_time()
    print('Total running time to load shape file: %s Seconds' % (end_time - start_time))
    
    return raw_network_csv


def linestring_to_points(feature,line):
    return {feature:line.coords}


def poly_to_points(feature,poly):
    return {feature:poly.exterior.coords}



def _loadNodes(network, raw_network):
    print('Loading nodes ...')
    _node_required_fields = {'A', 'B', 'geometry'}
    _node_optional_fields = {'zone_id'}
    
    fieldnames = list(raw_network)
    if '' in fieldnames:
        print('WARNING: columns with an empty header are detected in the network file. these columns will be skipped')
        fieldnames = [fieldname for fieldname in fieldnames if fieldname]
    
    fieldnames_set = set(fieldnames)
    for field in _node_required_fields:
        if field not in fieldnames_set:
            sys.exit(f'ERROR: required field ({field}) for generating node file does not exist in the network file')
            
    other_fields = list(fieldnames_set.difference(_node_required_fields.union(_node_optional_fields)))
    
    node_id_list = []
    node_coord_list = []
    node_dict = {}
    from_to_node_field = ['A', 'B']
    for index in raw_network.index:
        
        for i in range(2):
            node = Node(int(raw_network[from_to_node_field[i]][index]))
            
            if node.node_id in node_id_list:
                continue

            node_coords = raw_network["geometry"][index].coords
            node.geometry = geometry.Point(node_coords[i])
            node.x_coord, node.y_coord = float(node_coords[i][0]), float(node_coords[i][1])
            
            if node.node_id < 10000:
                node.zone_id = node.node_id
                
            # others
            for field in other_fields:
                node.other_attrs[field] = raw_network[field][index]
            
            node_dict[node.node_id] = node
            node_id_list.append(node.node_id)

    network.node_dict = node_dict
    network.node_other_attrs = other_fields
    
    print('%s nodes loaded' % len(node_id_list))
         
 

def _loadLinks(network, raw_network):
    print('Loading links ...')
    
    _link_required_fields = {'ID', 'A', 'B', 'DISTANCE', 'AMLANE', 'MDLANE', 'PMLANE', 'NTLANE', 'geometry'}
    _link_optional_fields = {}
    
    fieldnames = list(raw_network)
    if '' in fieldnames:
        print('WARNING: columns with an empty header are detected in the network file. these columns will be skipped')
        fieldnames = [fieldname for fieldname in fieldnames if fieldname]
    
    fieldnames_set = set(fieldnames)
    for field in _link_required_fields:
        if field not in fieldnames_set:
            sys.exit(f'ERROR: required field ({field}) for generating link file does not exist in the network file')
            
    other_fields = list(fieldnames_set.difference(_link_required_fields.union(_link_optional_fields)))    
    
    node_dict = network.node_dict
    link_dict = {}
    for index in raw_network.index:
        link = Link(int(index + 1))
        link.org_link_id = int(raw_network['ID'][index])
        
        from_node_id, to_node_id = int(raw_network['A'][index]), int(raw_network['B'][index])
        if from_node_id == to_node_id:
            print(f'WARNING: from_node and to_node of link {link.link_id} are the same')
            continue

        try:
            link.from_node = node_dict[from_node_id]
        except KeyError:
            print(f'WARNING: from_node {from_node_id} of link {link.link_id} does not exist in the node file')
            continue
        try:
            link.to_node = node_dict[to_node_id]
        except KeyError:
            print(f'WARNING: to_node {to_node_id} of link {link.link_id} does not exist in the node file')
            continue
            
        
        alpha_dict = {0:0.1,1:0.87,2:0.96,3:0.96,4:0.96,5:0.87,6:0.96} # classified by facility type
        beta_dict = {0:2,1:5,2:2.3,3:2.3,4:2.3,5:5,6:2.3} # classified by facility type
   
        
        time_period_dict = {1:'AM',2:'MD',3:'PM',4:'NT'}
        toll_allowed_uses_dict= {
            0: ['VDF_tollsov','VDF_tollhov2', 'VDF_tollhov3', 'VDF_tolltrk', 'VDF_tollapv', 'VDF_tollcom'],
            1: ['VDF_tollsov','VDF_tollhov2', 'VDF_tollhov3', 'VDF_tolltrk', 'VDF_tollapv', 'VDF_tollcom'],
            2: ['VDF_tollhov2', 'VDF_tollhov3'],
            3: ['VDF_tollhov3'],
            4: ['VDF_tollsov','VDF_tollhov2', 'VDF_tollhov3', 'VDF_tollapv', 'VDF_tollcom'],
            5: ['VDF_tollapv']
        }
        
        allowed_uses_dict = {0:'sov;hov2;hov3;trk;apv;com',1:'sov;hov2;hov3;trk;apv;com',2:'hov2;hov3',3:'hov3',
                         4:'sov;hov2;hov3;com;apv',5:'apv',6:'',7:'',8:'',9:'closed'}
        
        speed_class_dict = {0:17,1:17,2:17,3:23,4:29,5:35,6:40,11:63,12:63,13:69,14:69,15:75,16:75,21:40,	
                           22:40,23:52,24:52,25:58,26:58,31:40,32:40,33:46,34:46,35:46,36:52,41:35,42:35,43:35,	
                           44:40,45:40,46:40,51:52,52:52,53:58,54:58,55:58,56:63,61:23,62:23,63:35,64:35,65:40,66:58}
        
        capacity_class_dict = {0:3150,1:3150,2:3150,3:3150,4:3150,5:3150,6:3150,
                               11:1900,12:1900,13:2000,14:2000,15:2000,16:2000,21:600,22:800,23:960,24:960,
                               25:1100,26:1100,31:500,32:600,33:700,34:840,35:900,36:900,41:500,42:500,
                               43:600,44:800,45:800,46:800,51:1100,52:1200,53:1200,54:1400,55:1600,56:1600,
                               61:1000,62:1000,63:1000,64:1000,65:2000,66:2000}
        
        PHF_dict = {'AM':2.39776486,'MD':5.649424854,'PM':3.401127052,'NT':6.66626961}
            
        length_in_mile = raw_network['DISTANCE'][index]
        length = length_in_mile * 1609.34
        link.other_attrs['length_in_mile'] = float(length_in_mile) if length_in_mile else 0
        link.length = float(length) if length else 0
        
        
        lanes = raw_network['AMLANE'][index]
        link.lanes = int(lanes) if lanes else 0
        link.geometry = raw_network['geometry'][index]
        link_type = raw_network['ATYPE'][index] * 100 + raw_network['FTYPE'][index]
        vdf_code = link_type
        if link_type: link.link_type = int(link_type)
        if vdf_code: link.vdf_code = int(vdf_code)
        
        cap_class = int(raw_network['CAPCLASS'][index])
        capacity = capacity_class_dict[cap_class]
        link.capacity = int(capacity) if capacity else 0
            
        spd_class = int(raw_network['SPDCLASS'][index])
        free_speed = speed_class_dict[spd_class]
        link.free_speed = int(free_speed) if free_speed else 0
        
        FT = raw_network['FTYPE'][index]
        link.other_attrs['FT'] = FT
        
        AT = raw_network['ATYPE'][index]
        link.other_attrs['AT'] = FT
        
        for t_seq, t_period in time_period_dict.items():
            
            link.other_attrs['lanes'+str(t_seq)] = int(raw_network[str(t_period) + 'LANE' ][index])
            
            link.other_attrs['free_speed'+str(t_seq)] = link.free_speed
            link.other_attrs['VDF_PHF'+str(t_seq)] = float(PHF_dict[t_period])
            
            toll = raw_network[t_period + 'TOLL'][index]/100 # cents -> dollars
            allowed_uses_name = str('VDF_allowed_uses'+str(t_seq))
            
            try:
                allowed_uses_key = int(raw_network[str(t_period + 'LIMIT')][index])
            except ValueError:
                allowed_uses_key = 0
            
            allowed_uses = allowed_uses_dict[allowed_uses_key]
            
            for allowed_agent in toll_allowed_uses_dict[0]:
                link.other_attrs[str(allowed_agent) + str(t_seq)] = 0

            if allowed_uses_key >= 0: 
                link.other_attrs[allowed_uses_name] = allowed_uses
                if allowed_uses_key < 6:
                    for allowed_agent in toll_allowed_uses_dict[allowed_uses_key]:
                        link.other_attrs[str(allowed_agent) + str(t_seq)] = float(toll)
                
                     
            VDF_fftt = 60 * link.other_attrs['length_in_mile'] / link.free_speed
            if VDF_fftt: link.other_attrs['VDF_fftt'+str(t_seq)] = float(VDF_fftt)
            
            VDF_cap = link.other_attrs['lanes'+str(t_seq)] * link.capacity * link.other_attrs['VDF_PHF'+str(t_seq)]
            link.other_attrs['VDF_cap'+str(t_seq)] = float(VDF_cap) if VDF_cap else 0
            
            VDF_plf = 1/(float(PHF_dict[t_period]))
            link.other_attrs['VDF_plf'+str(t_seq)] = float(VDF_plf) if VDF_plf else 0
            
            
            
            VDF_alpha = float(alpha_dict[FT])
            link.other_attrs['VDF_alpha'+str(t_seq)] = float(VDF_alpha) if VDF_plf else 0
            
            VDF_beta = float(beta_dict[FT])
            link.other_attrs['VDF_beta'+str(t_seq)] = float(VDF_beta) if VDF_plf else 0
            
                
        for field in other_fields:
            link.other_attrs[field] = raw_network[field][index]
            

        link_dict[link.link_id] = link
        
    network.link_dict = link_dict
    network.link_other_attrs = other_fields
    
    print('%s links loaded' % len(link_dict))
         
         

def _buildnet(shapfile_path):
    
    raw_network = loadCSVFromSHP(shapfile_path)

    network = Network()
    _loadNodes(network, raw_network)
    _loadLinks(network, raw_network)

    return network



def _outputNode(network, output_folder):
    print('Generating node file ...')
    
    node_filename = 'node.csv'
    node_filepath = os.path.join(output_folder, node_filename)
    
    outfile = open(node_filepath, 'w', newline='',errors='ignore')

    writer = csv.writer(outfile)
    
    writer.writerow(['node_id', 'x_coord', 'y_coord', 'centroid', 'zone_id', 'geometry'])
    for node_id,node in network.node_dict.items():
        line = [node.node_id, node.x_coord, node.y_coord, node.centroid, node.zone_id, node.geometry]
        
        writer.writerow(line)
    outfile.close()
    print('node.csv generated')
    

    
def _outputLink(network, output_folder):
    print('Generating link file ...')
        
    link_filename = 'link.csv'
    link_filepath = os.path.join(output_folder, 'link.csv')
    
    outfile = open(link_filepath, 'w', newline='',errors='ignore')
        
    writer = csv.writer(outfile)
    link_header = ['first', 'link_id', 'from_node_id', 'to_node_id', 'lanes', 'length', 'dir_flag', 'geometry', 
                   'free_speed', 'capacity', 'link_type', 'vdf_code']
    
    first_link = network.link_dict[1]
    other_link_header = list(first_link.other_attrs.keys())
    
    link_header.extend(other_link_header)
    
    
    writer.writerow(link_header)
    
    for link_id, link in network.link_dict.items():
        
        line = ['', link.link_id, link.from_node.node_id, link.to_node.node_id, link.lanes, link.length, 
                link.dir_flag, link.geometry, link.free_speed, link.capacity, link.link_type, link.vdf_code]
        
        other_link_att_values = list(link.other_attrs.values())
        
        line.extend(other_link_att_values)
        
        writer.writerow(line)
    outfile.close()
    print('link.csv generated')
    






    
if __name__ == '__main__':
    
    shapfile_path = r'C:\Users\...\cube_net'
    output_folder = shapfile_path
    
    network = _buildnet(shapfile_path)
    _outputNode(network, output_folder)
    _outputLink(network, output_folder)
        

