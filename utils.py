# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 20:40:39 2021

@author: Peter

Function: This is a python file that contains most of the simple helpful
functions in our code.
"""

import os
import json


def extractFolderList(ab_path_now, data_folder):
    """
    Load the folder's path with various kinds of data
    Varibales:
        ab_path_now: The absolute path of current folder
        data_folder: The specified dataset:
            diffusiveWall_Jinkai: This is a dataset containing all txt files
                that give basic empirical data for completely diffuse
                scattering kernels.
            Maxwellian_Jinkai: This is a dataset for Maxwell-type scattering
                kernels.
    Return:
        folder_path_list: This is a list containing all the pathes of folders
            in the parent folder.
        name_list: This is a list giving the name of the child folders, but
            here we just use the velocity type without the postfix in the
            original folder name, thus you need to be careful to use this
            variable.
    """
    folder_path_list = []
    name_list = []

    data_path = os.path.join(ab_path_now, data_folder)
    path_dir = os.listdir(data_path)
    for child_dir in path_dir:
        child_path = os.path.join(data_path, child_dir)
        if os.path.isdir(child_path):
            # Obtain the folder path
            folder_path_list.append(child_path)

            # Extract the folder name
            name, _ = child_dir.split('_')
            name_list.append(name)
        elif os.path.isfile(child_path):
            print(child_path, 'is a file')
        else:
            print('It is a sepcial file')
    return folder_path_list, name_list


def extractFileInfo(folder_path):
    """
    Create a dictionary containing all the data, the key is the name of file,
    the value is the path of the file.
    Varibales:
        folder_path: This is a first child folder path, which chould be
            obtained from folder_path_list in extractFolderList function.
    Returns:
        file_name_list: All file names in this parent folder.
        file_list: The list contains all paths of useful text files.
        data_info: A dictionary contains the information of the folder, the key
            is the name of text file, the value is the path of corresponding
            file.
    """
    file_name_list = []
    file_list = []
    data_info = {}

    path_dir = os.listdir(folder_path)
    for child_dir in path_dir:
        name, _ = child_dir.split('.')
        if name and (name != 'REDEME'):
            file_name_list.append(name)
            file_path = os.path.join(folder_path, child_dir)
            file_list.append(file_path)
            data_info[name] = file_path

    return file_name_list, file_list, data_info


def createDataDic(absolute_path, total_folder):
    """
    Give the file path in dataset.
    Variables:
        absolute_path: It is a string containing the abolute path of our total
            foler.
        total_folder: It is the name of our total dataset, herem, it is SK_Data.
    Return:
        data_dic: A dictionary contatins all the information of the folder,
            the key is the name of text file, the value is the path of
            corresponding file.
    """
    # Get the absolute path of the current folder.
    # absolute_path = os.path.abspath('.')
    # print(absolute_path)
    # data_folder = 'SK_Data'  # This is the folder containing all the data we need
    folder_path_list, folder_name_list = extractFolderList(absolute_path,
                                                           total_folder)
    """
    print(
        "The folder path is {} \n, the folder name is {}.".format(folder_path_list,
                                                                  folder_name_list))
    """
    data_dic = {}  # data path dictionary
    for folder_path, folder_name in zip(folder_path_list, folder_name_list):
        _, _, data_info = extractFileInfo(folder_path)
        data_dic[folder_name] = data_info

    # store all information of our data into a json file
    json_str = json.dumps(data_dic, indent=4)
    with open('data_info.json', 'w') as json_file:
        json_file.write(json_str)

    return data_dic


if __name__ == '__main__':
    absolute_path = os.path.abspath('.')
    data_folder = 'SK_Data'
    data_dic = createDataDic(absolute_path, data_folder)
    print('-'*10, 'MAIN', '-'*10)
    print(data_dic.keys())