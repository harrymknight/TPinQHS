FROM ubuntu:22.04
RUN apt update && apt install python3-pip -y
COPY ./ ./
RUN pip3 install ./
CMD ["python3", "./app.py"]