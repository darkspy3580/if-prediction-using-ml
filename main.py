import streamlit as st
import base64

# Function to convert image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as img_file:
        img_base64 = base64.b64encode(img_file.read()).decode('utf-8')
    return img_base64

# Load the local image and convert it to base64
image_path = 'assets/home_bg.png'  # Replace with your local PNG image path
image_base64 = image_to_base64(image_path)
