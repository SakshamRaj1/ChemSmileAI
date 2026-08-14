import os
import time
import threading
from flask import Flask, jsonify

app = Flask(__name__)

def delayed_delete(file_path):
    # Wait for 10 minutes (10 minutes * 60 seconds = 600 seconds)
    time.sleep(600)
    if os.path.exists(file_path):
        os.remove(file_path)

@app.route('/create-json')
def create_json():
    file_path = 'data.json'
    
    # Write your JSON file logic here
    with open(file_path, 'w') as f:
        f.write('{"status": "temporary"}')

    # Start background thread to delete the file after 10 mins
    deletion_thread = threading.Thread(target=delayed_delete, args=(file_path,))
    deletion_thread.daemon = True
    deletion_thread.start()

    return jsonify({"message": "File created and will be deleted in 10 minutes."})

if __name__ == '__main__':
    app.run(debug=True)
