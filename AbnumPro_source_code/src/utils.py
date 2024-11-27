import xml.etree.ElementTree as ET

# 读取SVG文件
tree = ET.parse('./assets/images/imgt_number/imgt_num.svg')
import os
data = {
    "E": "1-1",
    "V": "1-2",
    "Q": "1-3",
    "L": "1-4",
    "S": "1-5",
    "G": "1-6",
    "A": "1-7",
    "R": "1-8",
    "K": "1-9",
    "T": "1-10",
    "I": "1-11",
    "H": "1-12",
    "P": "1-13",
    "W": "1-14",
    "Y": "1-15",
    "F": "1-16",
    "M": "1-17",
    "D": "1-18",
    "N": "1-19",
    "C": "2-1",
    "-" : "1-0 "
}
isolate_list = [4, 12, 13, 19, 21, 23, 35, 39, 41, 50, 52, 53, 55, 56, 71, 76, 78, 80, 87, 89, 91, 94, 100, 101, 104, 105, 106, 118]                            
aa_list = ["G","A","V","L","I","M","F","W"]
def clear_folder(folder_path):
    # 遍历文件夹中的所有文件和子文件夹
    for file in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file)
        # 如果是文件，则删除
        if os.path.isfile(file_path):
            os.remove(file_path)
        # 如果是文件夹，则递归清空
        elif os.path.isdir(file_path):
            clear_folder(file_path)

def getsvg(results):
    clear_folder('./assets/images/imgt_number/result')
    numbering, alignment_details, hit_tables = results[0]
    # num_list = numbering[0][0][0]
    # a = numbering[1][0][0]
    name_list = results[1][1]
    for i in range(0, len(numbering)):
        write_svg(name_list[i],numbering[i][0][0],str(i))
    
def write_svg(name,number_result,index):
    

    # 解析 SVG 文件
    tree = ET.parse('./assets/images/imgt_number/example_128.svg')
    if len(number_result)==127:
        tree = ET.parse('./assets/images/imgt_number/example_127.svg')
    root = tree.getroot()
    for i in range(0,len(number_result)):
        id_n = "n"+str(number_result[i][0][0])
        id_c = "c"+str(number_result[i][0][0])
        href_value = "#glyph"+data[number_result[i][1]]
       
        # element_with_id.attrib['ns1:href'] = href_value
        if number_result[i][1]=="-":
            if int(number_result[i][0][0])>26 and int(number_result[i][0][0])<39:
                 href_value = "#glyph"+"4-1"
            elif int(number_result[i][0][0])>55 and int(number_result[i][0][0])<66:
                href_value = "#glyph"+"4-2"
            elif int(number_result[i][0][0])>104 and int(number_result[i][0][0])<118:
                href_value = "#glyph"+"4-3"
            else:
                href_value = "#glyph"+"4-0"
        element_with_id = root.find(".//*[@id='"+id_c+"']")
        if number_result[i][1]=="P":
            
            element_with_id.set('style', "stroke:none;fill-rule:nonzero;fill:rgb(100%,100%,59.959412%);fill-opacity:1;")
        elif number_result[i][0][0] in isolate_list and number_result[i][1] in aa_list:
            element_with_id.set('style', "stroke:none;fill-rule:nonzero;fill:rgb(65.429688%,85.546875%,100%);fill-opacity:1;")
        else:
            element_with_id.set('style', " stroke:none;fill-rule:nonzero;fill:rgb(100%,100%,100%);fill-opacity:1;")

        element_with_id = root.find(".//*[@id='"+id_n+"']")
        element_with_id.set('{http://www.w3.org/1999/xlink}href', href_value)
    # 通过 ID 查找元素并修改属性
    # element_with_id = root.find(".//*[@id='n1']")
    # element_with_id.attrib['fill'] = 'red'  # 修改属性

    # 保存修改后的 SVG 文件
    tree.write("./assets/images/imgt_number/result/"+index+'.svg')

