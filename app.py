from web_app import app

# Point d'entrée pour Vercel serverless - ne pas modifier
if __name__ == "__main__":
    app.run(debug=True)
