## Hosting locally
To host the docs locally, navigate to the source directory and run
```
uv sync --group mkdocs
```

Then, run the following command to start the server:
```        # Go to the root directory
mkdocs serve
```

## Deploy
To deploy the docs, run
```
uv run mkdocs build --strict
pnpm dlx wrangler@latest deploy
```
