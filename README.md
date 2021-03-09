# ePotku

ePotku is a web version of [Potku](https://github.com/JYU-IBA/potku). It is used for time-of-flight elastic recoil detection analysis (ToF-ERDA).

The server can be run on Windows, Linux and macOS. Supported browsers are Firefox and Chrome.

# Getting started

Instructions for installing the backend and frontend are listed in their respective subfolders' READMEs. Starting with backend instructions is recommended.

# About Project ePotku
ePotku project developed a single-page web application (SPA) version of [Potku](https://github.com/JYU-IBA/potku) at the University of Jyväskylä in 2020. Potku (Finnish: kick, recoil), originally developed as a software project in 2013, is an open source analysis software for time-of-flight elastic recoil detection analysis (ToF-ERDA) measurement data. In ToF-ERDA, focused ion beams are used to investigate samples' elemental compositions. It was further expanded with a Monte-Carlo ERD simulation functionality by the Monisiro project in 2018.

ePotku can be used to create depth profiles, energy spectra and composition changes from ToF-ERDA measurement data. Computationally intensive calculations are offloaded to the server, speeding up analysis.

ePotku is also portable: one can work on a single request on multiple devices without having to manually move files. This also enables easy sharing of open data.

The front-end side ePotku was created with JavaScript, Vue.js and Plotly.js. The back-end utilizes Python 3.7, Flask, SQL and existing Potku code.

Project group members
Minja Hänninen
Ilari Jalli
Ville Kuokkanen
Pasi Niininen
Tuomas Pitkänen

Customer representatives
Mikko Laitinen
Timo Sajavaara
Kai Arstila
Jaakko Julin

Course instructors
Jonne Itkonen
Juhani Sundell
