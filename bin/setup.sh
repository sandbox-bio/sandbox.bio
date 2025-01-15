#!/bin/bash
# Setup script to contribute new tutorials to sandbox.bio

# Check that user has the right tools installed
echo -n "Checking for npm and node... "
if ! npm -v || ! node -v; then
    echo "Error: npm or node is not installed on your machine. Please follow instructions at https://docs.npmjs.com/downloading-and-installing-node-js-and-npm"
    exit
fi

# Set up mock .env file
echo "Setting up env file..."
cat > .env <<EOF
PUBLIC_USE_PRD_ASSETS=true
PUBLIC_SUPABASE_URL=https://127.0.0.1
PUBLIC_SUPABASE_API_KEY=mock
SUPABASE_API_KEY=mock
EOF

# Build data/ folder
echo "Building static assets..."
if ! ./bin/build.sh; then
    echo "Building assets failed."
    exit
fi

# Install dependencies
echo "Installing npm dependencies..."
if ! npm install > /dev/null; then
    echo "npm install failed."
    exit
fi

# Launch web server
npm run dev
