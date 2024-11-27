import os
import sys
from dataclasses import dataclass # pip install dataclasses
import pathlib
import jinja2
import json
from flask import Flask, render_template, request, send_from_directory
import webview # pip install pywebview
import pyglet
from src.api import API

if getattr(sys, "frozen", False):
    BASE_DIR = pathlib.Path(os.path.dirname(sys.executable))
elif __file__:
    BASE_DIR = pathlib.Path(os.path.dirname(__file__))

app = Flask(__name__, static_folder=str(BASE_DIR / 'assets'), template_folder=str(BASE_DIR / 'templates'))

@app.route("/")
def home():
    return render_template("index.html", context={"name": "PyDesktop App"}) 
@app.route("/index.html")
def move():
    return render_template("index.html", context={"name": "PyDesktop App"}) 
@app.route("/pages/forms/basic-forms.html")
def baic():
    return render_template("/pages/forms/basic-forms.html", context={"name": "PyDesktop App"}) 

@app.route("/switch_page", methods=['GET'])
def switch_page():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return json.dumps({'success': True, 'redirect_url': '/another_page'}), 200, {'ContentType':'application/json'}
@app.route("/number.html", methods=['GET'])
def switch_number():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return render_template("number.html")
    # return json.dumps({'success': True, 'redirect_url': '/'}), 200, {'ContentType':'application/json'}

@app.route("/number_result.html", methods=['GET'])
def number_result():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return render_template("number_result.html")
    # return json.dumps({'success': True, 'redirect_url': '/'}), 200, {'ContentType':'application/json'}

@app.route("/cdr.html", methods=['GET'])
def switch_cdr():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return render_template("cdr.html")
@app.route("/cdr_result.html", methods=['GET'])
def switch_cdr_result():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return render_template("cdr_result.html")
@app.route("/abr.html", methods=['GET'])
def switch_abr():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return render_template("abr.html")
@app.route("/abr_result.html", methods=['GET'])
def switch_abr_result():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return render_template("abr_result.html")

@app.route("/switch_page1", methods=['GET'])
def switch_page1():
    # 这里可以进行一些其他的操作或者逻辑

    # 返回一个JSON响应，包含重定向到另一个页面的URL
    return json.dumps({'success': True, 'redirect_url': '/'}), 200, {'ContentType':'application/json'}

@app.route("/another_page")
def another_page():
    return render_template("left_navigate.html")  # 渲染另一个页面
@app.route("/static/")
def serve_static():
    path_arg = request.args.get('path') # ?path=
    if path_arg:
        abs_path = pathlib.Path(path_arg).resolve()
        if abs_path.exists():
            file_dir = abs_path.parent
            file_name = abs_path.name
            return send_from_directory(file_dir, file_name)
    return render_template("index.html", context={"name": "PyDesktop App"}) 





if __name__ == "__main__":
    js_api = API(name='Justin')
    window_args = {
        "js_api": js_api,
        "width": 1800
    }

    
    window = webview.create_window("抗体的抗原结合区域计算分析", app,**window_args)
    js_api._window = window
    webview.start(debug=True, private_mode=False)