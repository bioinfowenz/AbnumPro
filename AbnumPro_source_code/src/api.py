from dataclasses import dataclass # pip install dataclasses
import pathlib
import json

import webview # pip install pywebview
from src.utils import getsvg
from num.abn import abn
from num.abr import abr
import csv
@dataclass
class API:
    name: str
    _window = None
    residue = ["G","A","V","L","I","P","F","Y","W","S","T","C","M","N","Q","D","E","K","R","H"]

    def selectImageFile(self):
        file_types = ("Images (*.png;*.jpg)", )
        return self._window.create_file_dialog(webview.OPEN_DIALOG, file_types=file_types)

    def selectFolder(self):
        return self._window.create_file_dialog(webview.FOLDER_DIALOG)

    def saveFile(self):
        file_types = ("Text (*.txt)", )
        return self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename='abc.txt')
    
    def saveProjectFile(self):
        file_types = ("PyDestkop (*.pydesktop)", )
        return self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename='abc.pydesktop')

    def triggerCreateThumbs(self, sourceDir):
        print('sourceDir', sourceDir)
        source_dir = sourceDir
        if isinstance(source_dir, list):
            source_dir = sourceDir[0]
        output_dir = dir_to_thumbnails(source_dir)
        return output_dir

    def myAPIRequest(self, jsonData):
        data  = json.loads(jsonData)
        destFilePath = data.get('filePath')
        print("my api request", jsonData, )
        if destFilePath:
            path = pathlib.Path(destFilePath).resolve()
            if path.exists():
                with open(path, 'w+') as f:
                    f.write(jsonData)


    def thisIsMyPyHandler(self, jsonData):
        print('this is my handler working', jsonData)

    def defaultHandleForm(self, jsonData):
        print("raw json data", jsonData)
        data = {}
        try:
            data = json.loads(jsonData)
        except:
            pass
        print("python data", data)
        print('myinputvalue', data.get('myinputname'))

    def triggerSomeJS(self):
        print("Trigger is working")
        print("run some long process")
        self._window.evaluate_js("helloWorldFromPy()")
        return 


    def getNumber(self, sequence,schema):
        print('sequence', sequence)
        print('schema', schema)
        numbered = "number_result"
        results = self.compute_number(sequence,schema)
        # source_dir = sourceDir
        # if isinstance(source_dir, list):
        #     source_dir = sourceDir[0]
        # output_dir = dir_to_thumbnails(source_dir)
        if schema=="Imgt":
            getsvg(results)
        return results
    def getNumberFile(self,sequencePath,schema):
        scheme = str.lower(schema)
        seq_dict = self.read_fasta(sequencePath)
        listseq = []
        listname = []
        sequeces = []
        for k,v in seq_dict.items():
            sequeces.append((k,v))
            listseq.append(v)
            listname.append(k)
        results = abn(sequeces, scheme, output=False)
        # print("sequence",sequenceFile)
        # results = self.getNumber(sequenceFile,schema)
        # print("sequence",sequenceFile)
        seq = [listseq, listname]
        return results,seq
    def selectFastaFile(self):
        file_types = ("Images (*.fasta;*.fa)", )
        return self._window.create_file_dialog(webview.OPEN_DIALOG, file_types=file_types)
    def compute_number(self,seq,schema):
        
        scheme = str.lower(schema)
        seqlin = seq.split("\n")
        listseq = []
        listname = []
        i=0
        seqnum = 0
        while(i<len(seqlin)):
            if(seqlin[i].startswith(">")):
                if(i+1<len(seqlin) and self.isantibody(seqlin[i+1])):
                    listname.append(seqlin[i])
                    i+=1
                    listseq.append(seqlin[i])
            elif(self.isantibody(seqlin[i])):
                listname.append("seq"+str(seqnum))
                listseq.append(seqlin[i])
                seqnum+=1

            i+=1
        sequeces = []
        for i in range(0,len(listseq)):
            
            sequeces.append((listname[i],listseq[i]))
        #
        #
        # for s in seqlin:
        #     if (s.startswith(">")):
        #         listname.append(s)
        #     else:
        #         listseq.append(s)
        seq = [listseq, listname]

        results = abn(sequeces, scheme, output=False)

        return results,seq

    def read_fasta(self,file_path):
        sequences = {}
        current_sequence = None

        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()  # 去除行尾的换行符和空格

                if line.startswith('>'):
                    # 如果行以 '>' 开头，表示新的序列开始
                    current_sequence = line[1:]  # 去除 '>' 字符
                    sequences[current_sequence] = ''
                elif current_sequence is not None:
                    # 如果不是 '>' 行，且当前有一个序列
                    sequences[current_sequence] += line

        return sequences    
    def isantibody(self,sequence):
        if(len(sequence)==0):
            return False
        for i in range(len(sequence)):
            if sequence[i] not in self.residue:
                return False
        return True
    def getCdrMessage(self, seq,scheme):
        # print("Received Message:", message)
        # seqlist = message["seq"]
        # scheme = message['scheme']
        scheme = str.lower(scheme)
        # sequences = []
        # for i in range(0, len(seqlist[0])):
        #     sequences.append((seqlist[1][i], seqlist[0][i]))
        # Hand the list of sequences to the abn function. Number them with the IMGT scheme
        # results = abn(seqlist, scheme="imgt", output=False)
        scheme1 = scheme
        if(scheme=="contact"):
            scheme1 = "martin"
        # results = abn(sequences, scheme1, output=False)
        results,seqArr = self.compute_number(seq,scheme1)
        # Unpack the results. We get three lists
        numbering, alignment_details, hit_tables = results

        cdr_list = []
        for i in range(0, len(numbering)):
            cdr_list.append(self.getCdr(scheme, numbering[i][0][0], alignment_details[i]))

        print(numbering)
        print(cdr_list)
        return seqArr,scheme,alignment_details,cdr_list

        # self.vBoxLayout.removeWidget(self.num)
        # # self.num = None
        # self.num.clearvis()
        # self.showCdr(cdr_list,numbering,alignment_details,scheme)
    def getCdrFile(self,sequencePath,schema):
        scheme = str.lower(schema)
        seq_dict = self.read_fasta(sequencePath)
        listseq = []
        listname = []
        sequeces = []
        for k,v in seq_dict.items():
            sequeces.append((k,v))
            listseq.append(v)
            listname.append(k)
        results = abn(sequeces, scheme, output=False)
        # print("sequence",sequenceFile)
        # results = self.getNumber(sequenceFile,schema)
        # print("sequence",sequenceFile)
        numbering, alignment_details, hit_tables = results
        seqArr = [listseq, listname]
        cdr_list = []
        for i in range(0, len(numbering)):
            cdr_list.append(self.getCdr(scheme, numbering[i][0][0], alignment_details[i]))

        return seqArr,scheme,alignment_details,cdr_list
    def getCdr(self, scheme, numbering, alignment_detail):
        cdr_list = ["", "", "","","","",""]
        fr_list = ["","","",""]
        heavy_index_range = [31, 35, 50, 65, 95, 102]
        light_index_range = [24, 34, 50, 56, 89, 97]
        cdr_range = {
            "kabat": [heavy_index_range, light_index_range]
        }
        chothia_h = [26,32,52,56,95,102]
        chothia_L = [24,34,50,56,89,97]
        cdr_range['chothia'] = [chothia_h,chothia_L]
        imgt_l=[27,32,50,51,89,97]
        imgt_h = [26,33,51,56,93,102]
        cdr_range['imgt'] = [imgt_h,imgt_l]
        contact_l=[30,36,46,55,89,96]
        contact_h = [30,35,47,58,93,101]
        cdr_range["contact"] = [contact_h,contact_l]
        for i in range(0,6):
            if(i%2==0):
                cdr_range['kabat'][0][i]-=1
                cdr_range['chothia'][0][i]-=1
                cdr_range['imgt'][0][i]-=1
                cdr_range["contact"][0][i]-=1
                cdr_range['kabat'][1][i]-=1
                cdr_range['chothia'][1][i]-=1
                cdr_range['imgt'][1][i]-=1
                cdr_range["contact"][1][i]-=1
            # else:
            #     cdr_range['kabat'][0]-=1
            #     cdr_range['chothia'][0]
            #     cdr_range['imgt'][0]
            #     cdr_range["contact"][0]
            #     cdr_range['kabat'][0]
            #     cdr_range['chothia'][0]
            #     cdr_range['imgt'][0]
            #     cdr_range["contact"][0]


        index = 0
        index_f = 0
        type = 0
        if (alignment_detail[0]['chain_type'] == 'H'):
            type = 0
        elif (alignment_detail[0]['chain_type'] == 'L' or alignment_detail[0]['chain_type'] == 'K'):
            type = 1
        for i in range(0, len(numbering)):
            if(index==6):
                if numbering[i][1]=="-":
                    continue
                cdr_list[index] = cdr_list[index] + numbering[i][1]
            else:
                if (numbering[i][0][0] <= cdr_range[str.lower(scheme)][type][index]):
                    if numbering[i][1]=="-":
                        continue
                    cdr_list[index] = cdr_list[index] + numbering[i][1]
                else:
                    index+=1
                    if numbering[i][1]=="-":
                        continue
                    cdr_list[index] = cdr_list[index] + numbering[i][1]
        # for i in range(0, len(numbering)):
        #     if (numbering[i][0][0] >= cdr_range[str.lower(scheme)][type][index * 2]):
        #         if (numbering[i][0][0] <= cdr_range[str.lower(scheme)][type][index * 2 + 1]):
        #             cdr_list[index] = cdr_list[index] + numbering[i][1]
        #         else:
        #             index += 1
        #             if (index == 3):
        #                 break

        return cdr_list
    def saveNumberFile(self):
        file_types = ("Text (*.txt)", )
        return self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename='numberingResult.txt')
    def writeNumberResult(self,filePath,results):
        # scheme = self.SchemeCombox.currentText()
        # results = abn(sequences, scheme, output=False)
        numbering, alignment_details, hit_tables = results
        # filePath = self.outputPath
        file = open(filePath,"a")
        for i in range(0,len(numbering)):
            file.write(alignment_details[i][0]["query_name"]+"\n")
            str_num = ""
            str_ali = ""
            str_seq = ""
            for num,seq in numbering[i][0][0]:
                str1 = ""+str(num[0])+num[1]
                str2 = "|"
                str3 = ""+seq
                width = 5
                str_num = str_num+str1.center(width)
                str_ali = str_ali+str2.center(width)
                str_seq = str_seq+str3.center(width)
            file.write(str_num+"\n")
            file.write(str_ali+"\n")
            file.write(str_seq+"\n")
        file.close()
        return "保存成功"
    def saveCdrFile(self):
        file_types = ("Text (*.csv)", )
        return self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename='cdrResult.csv')
    def writeCdrResult(self,filePath,cdr_list,seqlist):
        # scheme = self.SchemeCombox.currentText()
        # results = abn(sequences, scheme, output=False)
        # numbering, alignment_details, hit_tables = results
        # filePath = self.outputPath
        # file = open(filePath,"a")
        # cdrs = []
        header = ["名称","FR-1","CDR-1","FR-2","CDR-2","FR-3","CDR-3","FR-4"]
        
        for i in range(0,len(cdr_list)):
            name = seqlist[1][i]
            # alignment_details[i][0]["query_name"]
            # cdr_list = self.getCdr(scheme,name,numbering[i][0][0], alignment_details[i])
            cdr_list[i].insert(0,name)
            # cdrs.append(cdr_list)
        cdr_list.insert(0,header)
        file = open(filePath, "a", newline="")
        writer = csv.writer(file)
        writer.writerows(cdr_list)
        file.close()
        return "保存成功"
    def getAbrMessage(self, seq):
        # print("Received Message:", message)
        # seqlist = message["seq"]
        # scheme = message['scheme']
        # scheme = str.lower(scheme)
        # sequences = []
        # for i in range(0, len(seqlist[0])):
        #     sequences.append((seqlist[1][i], seqlist[0][i]))
        # Hand the list of sequences to the abn function. Number them with the IMGT scheme
        # results = abn(seqlist, scheme="imgt", output=False)
        # scheme1 = scheme
        # if(scheme=="contact"):
        #     scheme1 = "martin"
        # results = abn(sequences, scheme1, output=False)
        sequences = self.getSeq(seq)
        listseq, listname = sequences
        seq = []
        for i in range(len(listseq)):
            seq.append((listname[i],listseq[i]))

        # Unpack the results. We get three lists
        results = abr(seq)

        # cdr_list = []
        # for i in range(0, len(numbering)):
        #     cdr_list.append(self.getCdr(scheme, numbering[i][0][0], alignment_details[i]))

        # print(numbering)
        # print(cdr_list)
        return sequences,results

        # self.vBoxLayout.removeWidget(self.num)
        # # self.num = None
        # self.num.clearvis()
        # self.showCdr(cdr_list,numbering,alignment_details,scheme)
    def getSeq(self,seq):
        seqlin = seq.split("\n")
        listseq = []
        listname = []
        i=0
        seqnum = 0
        while(i<len(seqlin)):
            if(seqlin[i].startswith(">")):
                if(i+1<len(seqlin) and self.isantibody(seqlin[i+1])):
                    listname.append(seqlin[i])
                    i+=1
                    listseq.append(seqlin[i])
            elif(self.isantibody(seqlin[i])):
                listname.append("seq"+str(seqnum))
                listseq.append(seqlin[i])
                seqnum+=1

            i+=1
        sequeces = []
        for i in range(0,len(listseq)):
            
            sequeces.append((listname[i],listseq[i]))
        #
        #
        # for s in seqlin:
        #     if (s.startswith(">")):
        #         listname.append(s)
        #     else:
        #         listseq.append(s)
        seq = [listseq, listname]
        return seq
    def getAbrFile(self,sequencePath):
        # scheme = str.lower(schema)
        seq_dict = self.read_fasta(sequencePath)
        listseq = []
        listname = []
        sequeces = []
        for k,v in seq_dict.items():
            sequeces.append((k,v))
            listseq.append(v)
            listname.append(k)
        results = abr(sequeces)
        # print("sequence",sequenceFile)
        # results = self.getNumber(sequenceFile,schema)
        # print("sequence",sequenceFile)
        # numbering, alignment_details, hit_tables = results
        seqArr = [listseq, listname]
        # cdr_list = []
        # for i in range(0, len(numbering)):
        #     cdr_list.append(self.getCdr(scheme, numbering[i][0][0], alignment_details[i]))

        return seqArr,results
    def saveAbrFile(self):
        file_types = ("Text (*.csv)", )
        return self._window.create_file_dialog(webview.SAVE_DIALOG, save_filename='abrResult.csv')
    def writeAbrResult(self,filePath,cdr_list,seqlist):
        # scheme = self.SchemeCombox.currentText()
        # results = abn(sequences, scheme, output=False)
        # numbering, alignment_details, hit_tables = results
        # filePath = self.outputPath
        # file = open(filePath,"a")
        # cdrs = []
        header = ["名称","FR-1","ABR-1","FR-2","ABR-2","FR-3","ABR-3","FR-4"]
        
        for i in range(0,len(cdr_list)):
            name = seqlist[1][i]
            # alignment_details[i][0]["query_name"]
            # cdr_list = self.getCdr(scheme,name,numbering[i][0][0], alignment_details[i])
            cdr_list[i].insert(0,name)
            # cdrs.append(cdr_list)
        cdr_list.insert(0,header)
        file = open(filePath, "a", newline="")
        writer = csv.writer(file)
        writer.writerows(cdr_list)
        file.close()
        return "保存成功"
