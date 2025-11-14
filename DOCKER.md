# Docker Deployment Guide

## Prerequisites

- Docker installed (version 20.10 or higher)
- Docker Compose installed (version 2.0 or higher)
- AWS CLI configured (for downloading data)
- At least 4GB of free disk space

## Quick Start

```bash
# 1. Clone the repository (if not already done)
git clone <repository-url>
cd rnaseq-portal

# 2. Download the data directory from S3
aws s3 cp s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/rnaseq/portal/rnaseq-portal-data.zip .
unzip rnaseq-portal-data.zip

# 3. Build and start the container
docker-compose up -d

# 4. View logs to confirm startup
docker-compose logs -f

# 5. Access the app
# Open http://localhost:3838 in your browser
```

## Architecture

### Dockerfile Details

The Docker image is built on `rocker/r-ver:4.3` which provides:
- R version 4.3
- Essential system libraries for R packages
- Optimized for running R applications

**Layers:**
1. Base R image
2. System dependencies (libcurl, libssl, libxml2, etc.)
3. CRAN packages (shiny, ggplot2, foreach, doMC, gridExtra)
4. Bioconductor packages (DESeq2)
5. Application code (main/ and scripts/)

**Data Strategy:**
- Data directory is mounted as a volume (not included in image)
- Keeps image size small (~1.5GB vs ~3GB with data)
- Allows easy data updates without rebuilding

### docker-compose.yml

Simplifies container management with:
- Port mapping: `3838:3838`
- Volume mounting: `./data:/app/data:ro` (read-only)
- Health checks: Ensures container is responding
- Restart policy: `unless-stopped` for high availability

## Common Commands

### Create certs
```bash
mkdir certs/; cd certs
openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout server.key -out server.crt
```

### Building

```bash
# Build or rebuild the image
docker-compose build

# Build with no cache (force fresh build)
docker-compose build --no-cache

# Pull latest base image and rebuild
docker-compose build --pull
```

### Running

```bash
# Start in detached mode (background)
docker-compose up -d

# Start with logs visible
docker-compose up

# Scale (not applicable for single service, but possible)
docker-compose up -d --scale rnaseq-portal=1
```

### Monitoring

```bash
# View logs
docker-compose logs

# Follow logs in real-time
docker-compose logs -f

# View last 100 lines
docker-compose logs --tail=100

# Check container status
docker-compose ps

# View resource usage
docker stats rnaseq-portal
```

### Stopping

```bash
# Stop container (preserves data)
docker-compose stop

# Stop and remove container
docker-compose down

# Stop, remove, and delete volumes (WARNING: destructive)
docker-compose down -v
```

### Maintenance

```bash
# Restart the container
docker-compose restart

# Execute commands inside running container
docker-compose exec rnaseq-portal R --version
docker-compose exec rnaseq-portal ls -la /app/data

# Get a shell inside the container
docker-compose exec rnaseq-portal bash
```

## Troubleshooting

### Image Build Fails

**Issue:** Package installation fails during build

```bash
# Check if you have a stable internet connection
ping cloud.r-project.org

# Try building with verbose output
docker-compose build --progress=plain

# Clear Docker cache and rebuild
docker system prune -a
docker-compose build --no-cache
```

### Container Exits Immediately

**Issue:** Container starts but immediately exits

```bash
# Check logs for error messages
docker-compose logs

# Common causes:
# 1. Data directory not mounted correctly
# 2. Missing required data files
# 3. R package errors

# Verify data directory exists and contains files
ls -la data/
ls -la data/gene-counts_FPE/ | head
```

### Port Already in Use

**Issue:** Error binding to port 3838

```bash
# Find what's using port 3838
lsof -i :3838  # macOS/Linux
netstat -ano | findstr :3838  # Windows

# Option 1: Stop the conflicting service
# Option 2: Use a different port
# Edit docker-compose.yml: change "3838:3838" to "3839:3838"
# Then access at http://localhost:3839
```

### App Loads Slowly

**Issue:** Initial page load takes 30-60 seconds

This is **normal behavior** - the app needs to:
1. Load count matrix (62,757 genes × 231 samples)
2. Load pre-computed DESeq results
3. Initialize all modules

**Solutions:**
- Wait for initial load (subsequent interactions are fast)
- Pre-compute more data (see README_PRECOMPUTE.md)
- Increase container resources in Docker Desktop

### Data Volume Not Mounting

**Issue:** App can't find data files

```bash
# Verify volume is mounted
docker inspect rnaseq-portal | grep -A 10 Mounts

# Check data directory path in docker-compose.yml
cat docker-compose.yml | grep -A 2 volumes

# Ensure data directory exists before starting
ls -la data/

# Try absolute path in docker-compose.yml
# Replace: ./data:/app/data:ro
# With: /absolute/path/to/data:/app/data:ro
```

## Production Deployment

### Using a Different Port

Edit `docker-compose.yml`:

```yaml
ports:
  - "8080:3838"  # Access at http://localhost:8080
```

### Setting Resource Limits

Edit `docker-compose.yml`:

```yaml
services:
  rnaseq-portal:
    # ... existing config ...
    deploy:
      resources:
        limits:
          cpus: '2.0'
          memory: 4G
        reservations:
          cpus: '1.0'
          memory: 2G
```

### Behind a Reverse Proxy (nginx)

Example nginx configuration:

```nginx
server {
    listen 80;
    server_name rnaseq.example.com;

    location / {
        proxy_pass http://localhost:3838;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 600s;
    }
}
```

### Environment Variables

Add to `docker-compose.yml`:

```yaml
environment:
  - SHINY_LOG_LEVEL=INFO  # DEBUG for verbose, WARN for less
  - SHINY_LOG_STDERR=1    # Log to stderr
  - R_MAX_VSIZE=8Gb       # Increase R memory limit
```

## Security Considerations

### Read-Only Data Volume

The data volume is mounted read-only (`:ro`) to prevent accidental modifications:

```yaml
volumes:
  - ./data:/app/data:ro  # read-only
```

### Running as Non-Root User

For enhanced security, run as non-root user. Edit `Dockerfile`:

```dockerfile
# Add before CMD
RUN useradd -m -u 1000 shiny
USER shiny
```

### Network Isolation

For production, use Docker networks:

```yaml
services:
  rnaseq-portal:
    networks:
      - rnaseq-net

networks:
  rnaseq-net:
    driver: bridge
```

## Performance Optimization

### Pre-loading at Build Time

If you want faster startup (but larger image), load data during build:

```dockerfile
# Add to Dockerfile before CMD
COPY data/ /app/data/
```

⚠️ **Warning:** This increases image size by ~2GB

### Multi-stage Build

For smaller production images, use multi-stage builds (advanced).

### Caching

Docker caches layers - order matters:
1. Install system packages (changes rarely)
2. Install R packages (changes rarely)
3. Copy application code (changes frequently)

## Health Checks

The container includes health checks. View status:

```bash
# Check health status
docker inspect --format='{{.State.Health.Status}}' rnaseq-portal

# View health check logs
docker inspect --format='{{range .State.Health.Log}}{{.Output}}{{end}}' rnaseq-portal
```

## Backup and Restore

### Backup Data

```bash
# Backup data directory
tar -czf rnaseq-data-backup-$(date +%Y%m%d).tar.gz data/
```

### Restore Data

```bash
# Extract backup
tar -xzf rnaseq-data-backup-20251104.tar.gz
```

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Build and Push Docker Image

on:
  push:
    branches: [main]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Build Docker image
        run: docker build -t rnaseq-portal:latest .
      - name: Test image
        run: |
          docker run -d --name test -p 3838:3838 rnaseq-portal:latest
          sleep 30
          curl -f http://localhost:3838 || exit 1
```

## Support

For issues specific to Docker deployment:
1. Check container logs: `docker-compose logs`
2. Verify data volume: `docker inspect rnaseq-portal`
3. Test R packages: `docker-compose exec rnaseq-portal R -e "library(DESeq2)"`
4. Check application logs inside container: `docker-compose exec rnaseq-portal cat /tmp/*.log`

---

**🐳 For general app documentation, see [README.md](README.md)**

