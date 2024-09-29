# t14repythonjupyter
## 環境構築


rye init t14re-python-jupyter
cd t14re-python-jupyter

rye pin 3.8
rye sync

rye add opencv-python
rye add scipy==1.7.1
rye add numpy==1.21.5
rye add pyserial==3.5
rye add jupyter==1.0.0
rye sync

# jupyter notebookの起動
jupyter notebook