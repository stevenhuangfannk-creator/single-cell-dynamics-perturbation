# 创建独立环境。不要把本脚本指向已有环境；Python 3.11 + torch 2.5.1 保留官方检查点兼容性。
param([string]$Environment = 'work/sctour/.venv', [switch]$CPU)
$ErrorActionPreference = 'Stop'
uv venv $Environment --python 3.11
$TaskPython = Join-Path $Environment 'Scripts/python.exe'
if ($CPU) { $TorchIndex = 'https://download.pytorch.org/whl/cpu' }
else { $TorchIndex = 'https://download.pytorch.org/whl/cu121' }
uv pip install --python $TaskPython torch==2.5.1 --index-url $TorchIndex
uv pip install --python $TaskPython 'torch==2.5.1' 'sympy==1.13.1' 'numpy==1.26.4' 'scipy==1.13.1' 'pandas==2.2.3' 'scanpy==1.9.8' 'anndata==0.10.8' 'matplotlib==3.8.4' 'scikit-learn==1.6.1' 'scikit-misc==0.5.3' 'torchdiffeq==0.2.5' 'sctour==1.0.0' 'psutil==7.2.2'
uv pip check --python $TaskPython
if ($LASTEXITCODE -ne 0) { throw 'scTour dependency compatibility check failed' }
& $TaskPython -c 'import torch,sctour,scanpy; print(torch.__version__, torch.cuda.is_available(), scanpy.__version__)'
