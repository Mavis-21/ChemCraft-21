import streamlit as st
import sqlite3
import requests
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

# --- SETUP & SESSION STATE ---
def initialisation():
    if "logged_in" not in st.session_state:
        st.session_state.logged_in = False
        st.session_state.guest = False
        st.session_state.user = None
        st.session_state.page = "home" # Default page
    if "guest_history" not in st.session_state:
        st.session_state.guest_history = []

# --- SQLITE DATABASE FUNCTIONS ---
def connection():
    # Creates/Connects to chemcraft.db
    con = sqlite3.connect('chemcraft.db', check_same_thread=False)
    return con, con.cursor()

def create_tables():
    con, cur = connection()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS users(
            userid INTEGER PRIMARY KEY AUTOINCREMENT,
            username TEXT UNIQUE,
            passwd TEXT,
            email TEXT,
            typ TEXT
        )""")
    con.commit()
    con.close()

def create_usertable(user):
    con, cur = connection()
    # Create a table for the user if it doesn't exist
    cur.execute(f'create table if not exists "{user}" (user text, searched text, smiles text)')
    con.commit()
    con.close()

def register_user(user, passwd, email, typ):
    try:
        con, cur = connection()
        cur.execute("INSERT INTO users (username, passwd, email, typ) VALUES (?, ?, ?, ?)", 
                    (user, passwd, email, typ))
        con.commit()
        con.close()
        return True
    except:
        return False

def check_login(user, passwd):
    con, cur = connection()
    try:
        cur.execute("SELECT passwd FROM users WHERE username = ?", (user,))
        res = cur.fetchone()
        if res and res[0] == passwd: return True
    except: pass
    finally: con.close()
    return False

def get_history(username):
    con, cur = connection()
    try:
        # Check if user table exists
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (username,))
        if not cur.fetchone(): return []
        
        cur.execute(f'SELECT searched FROM "{username}"')
        data = cur.fetchall()
        return [d[0] for d in data]
    except: return []
    finally: con.close()

def save_history(user, iupac, smiles):
    create_usertable(user)
    con, cur = connection()
    try:
        # Check duplicate
        cur.execute(f'SELECT * FROM "{user}" WHERE searched = ?', (iupac,))
        if not cur.fetchone():
            cur.execute(f'INSERT INTO "{user}" VALUES (?, ?, ?)', (user, iupac, smiles))
            con.commit()
    finally: con.close()

# --- CHEMISTRY FUNCTIONS ---
def iupac_to_smiles(name):
    try:
        url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        r = requests.get(url, timeout=5)
        if r.status_code == 200: return r.json().get("smiles")
    except: return None

def get_3d_block(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        return Chem.MolToMolBlock(mol)
    except: return None

# --- PAGES ---
def home_page():
    st.markdown("<h1 style='text-align: center; color: #4facfe;'>ChemCraft</h1>", unsafe_allow_html=True)
    st.markdown("<h4 style='text-align: center;'>Interactive 3D Molecule Visualizer</h4>", unsafe_allow_html=True)
    st.write("")
    
    c1, c2, c3 = st.columns(3)
    with c1:
        if st.button("Sign Up", use_container_width=True): st.session_state.page = "signup"
    with c2:
        if st.button("Log In", use_container_width=True): st.session_state.page = "login"
    with c3:
        if st.button("Guest Mode", use_container_width=True):
            st.session_state.guest = True
            st.session_state.page = "dashboard"

def signup_page():
    st.header("Sign Up")
    with st.form("signup"):
        u = st.text_input("Username")
        p = st.text_input("Password", type="password")
        e = st.text_input("Email")
        role = st.selectbox("Role", ["Student", "Teacher", "Other"])
        if st.form_submit_button("Create Account"):
            if register_user(u, p, e, role):
                st.success("Account created! Please log in.")
                st.session_state.page = "login"
                st.rerun()
            else:
                st.error("Username already taken.")
    if st.button("Back"): st.session_state.page = "home"

def login_page():
    st.header("Login")
    with st.form("login"):
        u = st.text_input("Username")
        p = st.text_input("Password", type="password")
        if st.form_submit_button("Log In"):
            if check_login(u, p):
                st.session_state.logged_in = True
                st.session_state.user = u
                st.session_state.page = "dashboard"
                st.rerun()
            else:
                st.error("Invalid credentials.")
    if st.button("Back"): st.session_state.page = "home"

def dashboard():
    # SIDEBAR NAVIGATION
    with st.sidebar:
        st.header(f"Hello, {st.session_state.user if st.session_state.user else 'Guest'}!")
        if st.button("🔍 New Search", use_container_width=True): 
            st.session_state.view_mode = "search"
        
        if st.session_state.logged_in:
            st.write("---")
            st.write("**History**")
            hist = get_history(st.session_state.user)
            for item in hist[-5:]:
                if st.button(f"🧪 {item}", key=item):
                    st.session_state.view_mode = "history"
                    st.session_state.current_mol = item
        
        st.write("---")
        if st.button("Log Out", use_container_width=True):
            st.session_state.logged_in = False
            st.session_state.guest = False
            st.session_state.page = "home"
            st.rerun()

    # MAIN CONTENT
    st.title("Molecule Viewer")
    
    # Determine what to show
    if "view_mode" not in st.session_state: st.session_state.view_mode = "search"
    
    target_mol = None
    
    if st.session_state.view_mode == "search":
        target_mol = st.text_input("Enter IUPAC Name (e.g., aspirin, benzene):")
    elif st.session_state.view_mode == "history":
        target_mol = st.session_state.get("current_mol")
        st.info(f"Viewing from history: {target_mol}")

    if target_mol:
        if st.button("Visualize") or st.session_state.view_mode == "history":
            with st.spinner("Generating 3D Model..."):
                smiles = iupac_to_smiles(target_mol)
                if smiles:
                    st.success(f"SMILES: {smiles}")
                    # Save to history if logged in
                    if st.session_state.logged_in:
                        save_history(st.session_state.user, target_mol, smiles)
                    
                    # Render 3D
                    blk = get_3d_block(smiles)
                    if blk:
                        view = py3Dmol.view(width=700, height=500)
                        view.addModel(blk, "mol")
                        view.setStyle({'stick': {}})
                        view.zoomTo()
                        st.components.v1.html(view._make_html(), height=500)
                    else:
                        st.error("Could not generate 3D structure.")
                else:
                    st.error("Molecule not found. Check spelling.")

# --- MAIN APP LOOP ---
create_tables()
initialisation()

if st.session_state.page == "home": home_page()
elif st.session_state.page == "signup": signup_page()
elif st.session_state.page == "login": login_page()
elif st.session_state.page == "dashboard": dashboard()
