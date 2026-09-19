FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# 先装依赖（利用层缓存）
COPY requirements.txt ./
RUN pip install -r requirements.txt

# 再拷代码
COPY app ./app
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh \
    && mkdir -p /data

# 数据库文件挂载点
VOLUME ["/data"]
ENV PRIMER_DB_PATH=/data/primer_dimer.db \
    PRIMER_HOST=0.0.0.0 \
    PRIMER_PORT=8000

EXPOSE 8000

# 依赖安装（构建期）与数据库初始化（启动期）都在容器流程内
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
