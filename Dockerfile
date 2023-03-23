FROM python:3.8.12
USER root

RUN python -m pip install numpy==1.20.3
RUN python -m pip install pandas==1.3.4
RUN python -m pip install beautifulsoup4==4.10.0
RUN python -m pip install requests==2.26.0