#!/bin/bash

# sphinx-apidoc -f -o ./source/ ../python_driver/
make html
chromium-browser build/html/index.html
